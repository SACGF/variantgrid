from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.models import Analysis
from analysis.models.enums import NodeMatchInput
from analysis.models.nodes.filters.classifications_node import ClassificationsNode
from analysis.models.nodes.node_counts import get_extra_filters_q
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import get_fake_annotation_version
from classification.enums import ClinicalSignificance, SomaticClinicalSignificance, SpecialEKeys, SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.grid_columns.custom_columns import get_variantgrid_extra_annotate
from snpdb.models import GenomeBuild, Variant
from snpdb.models.models_enums import AlleleOriginFilterDefault, BuiltInFilters
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


class ClassificationTestDataMixin:
    """ A published germline P plus somatic Tier I, Tier III and Tier I/II classifications """

    def setUp(self):
        super().setUp()
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.germline_p = self._create_classification("germline", clinical_significance="P")
        self.somatic_tier_1 = self._create_classification("somatic",
                                                          somatic_clinical_significance=SomaticClinicalSignificance.TIER_1)
        self.somatic_tier_3 = self._create_classification("somatic",
                                                          somatic_clinical_significance=SomaticClinicalSignificance.TIER_3)
        self.somatic_tier_1_or_2 = self._create_classification("somatic",
                                                               somatic_clinical_significance=SomaticClinicalSignificance.TIER_1_OR_2)

    def _create_classification(self, allele_origin: str, clinical_significance: str = None,
                               somatic_clinical_significance: str = None) -> Classification:
        data = {
            SpecialEKeys.C_HGVS: {'value': 'c.301A>C'},
            SpecialEKeys.ALLELE_ORIGIN: {'value': allele_origin},
        }
        if clinical_significance:
            data[SpecialEKeys.CLINICAL_SIGNIFICANCE] = {'value': clinical_significance}
        if somatic_clinical_significance:
            data[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE] = {'value': somatic_clinical_significance}
        classification = Classification.create(
            user=self.user,
            lab=self.lab,
            lab_record_id=None,
            data=data,
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)
        classification.publish_latest(user=self.user)
        return classification


@override_settings(CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class ClassificationsNodeFilterTest(ClassificationTestDataMixin, TestCase):
    """ The node's germline/somatic Q building - allele origin buckets paired with the right
        clinical significance field, all-ticked meaning no significance filter """

    def _matching_classification_ids(self, node: ClassificationsNode) -> set:
        cm_qs = ClassificationModification.latest_for_user(self.user, published=True)
        return set(cm_qs.filter(node._classification_match_q()).values_list("classification_id", flat=True))

    def test_default_matches_everything(self):
        node = ClassificationsNode()
        self.assertEqual(self._matching_classification_ids(node),
                         {self.germline_p.pk, self.somatic_tier_1.pk, self.somatic_tier_3.pk,
                          self.somatic_tier_1_or_2.pk})

    def test_germline_origin_excludes_somatic(self):
        node = ClassificationsNode(allele_origin=AlleleOriginFilterDefault.GERMLINE)
        self.assertEqual(self._matching_classification_ids(node), {self.germline_p.pk})

    def test_somatic_origin_excludes_germline(self):
        node = ClassificationsNode(allele_origin=AlleleOriginFilterDefault.SOMATIC)
        self.assertEqual(self._matching_classification_ids(node),
                         {self.somatic_tier_1.pk, self.somatic_tier_3.pk, self.somatic_tier_1_or_2.pk})

    def test_germline_significance_subset_keeps_somatic_side(self):
        node = ClassificationsNode(pathogenic=False)
        self.assertEqual(self._matching_classification_ids(node),
                         {self.somatic_tier_1.pk, self.somatic_tier_3.pk, self.somatic_tier_1_or_2.pk})

    def test_somatic_tier_subset_includes_tier_1_or_2_records(self):
        # There is no Tier I/II filter option - a record recorded as "Tier I/II" matches
        # whenever Tier I or Tier II is selected
        node = ClassificationsNode(allele_origin=AlleleOriginFilterDefault.SOMATIC,
                                   tier_2=False, tier_3=False, tier_4=False)
        self.assertEqual(self._matching_classification_ids(node),
                         {self.somatic_tier_1.pk, self.somatic_tier_1_or_2.pk})

    def test_somatic_tier_subset_without_tier_1_or_2(self):
        node = ClassificationsNode(allele_origin=AlleleOriginFilterDefault.SOMATIC,
                                   tier_1=False, tier_2=False, tier_4=False)
        self.assertEqual(self._matching_classification_ids(node), {self.somatic_tier_3.pk})

    def test_get_classifications_qs_ors_germline_and_somatic(self):
        # A record matches if either its germline or somatic clinical significance is in its list
        qs = Classification.get_classifications_qs(
            self.user,
            clinical_significance_list=[ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC],
            somatic_clinical_significance_list=SomaticClinicalSignificance.TIER_1_AND_2_VALUES)
        self.assertEqual(set(qs.values_list("pk", flat=True)),
                         {self.germline_p.pk, self.somatic_tier_1.pk, self.somatic_tier_1_or_2.pk})

    def test_chips_reflect_selection(self):
        node = ClassificationsNode(allele_origin=AlleleOriginFilterDefault.SOMATIC,
                                   tier_2=False, tier_3=False, tier_4=False)
        chip_texts = [chip.text for chip in node.get_node_chips()]
        self.assertEqual(chip_texts, ["Somatic", "I"])

    def test_default_node_has_no_chips(self):
        self.assertEqual(ClassificationsNode().get_node_chips(), [])


@override_settings(CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class InternalClassificationGridAnnotationTest(ClassificationTestDataMixin, TestCase):
    """ The grid box annotations - germline max, somatic max ordered by summary sort score
        (not lexically - tier_1 outranks tier_3) and the '|'-joined somatic summary """

    def setUp(self):
        super().setUp()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant = slowly_create_test_variant("1", 1000, "A", "T", grch37)
        allele = create_mock_allele(self.variant, grch37)
        for classification in [self.germline_p, self.somatic_tier_1, self.somatic_tier_3,
                               self.somatic_tier_1_or_2]:
            classification.allele = allele
            classification.save()

    def test_grid_annotations(self):
        annotations = get_variantgrid_extra_annotate(self.user)
        values = Variant.objects.filter(pk=self.variant.pk).annotate(**annotations).values(
            "max_internal_classification", "max_internal_somatic_classification",
            "internally_classified_somatic").get()
        self.assertEqual(values["max_internal_classification"], ClinicalSignificance.PATHOGENIC)
        self.assertEqual(values["max_internal_somatic_classification"], SomaticClinicalSignificance.TIER_1)
        self.assertEqual(set(values["internally_classified_somatic"].split("|")),
                         {"U", "tier_1", "tier_3", "tier_1_or_2"})


@override_settings(CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class BuiltInClassificationFilterTest(ClassificationTestDataMixin, TestCase):
    """ The classification built-in filters / node counts - Classified is origin-agnostic,
        Classified Pathogenic is germline LP/P, Classified Tier I/II is somatic """

    def setUp(self):
        super().setUp()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.annotation_version = get_fake_annotation_version(grch37)
        self.analysis = Analysis(genome_build=grch37)
        self.analysis.set_defaults_and_save(self.user)
        self.variant_pks = {}
        for position, classification in enumerate([self.germline_p, self.somatic_tier_1, self.somatic_tier_3,
                                                   self.somatic_tier_1_or_2], start=1000):
            variant = slowly_create_test_variant("1", position, "A", "T", grch37)
            classification.allele = create_mock_allele(variant, grch37)
            classification.save()
            self.variant_pks[classification.pk] = variant.pk

    def _filter_variant_pks(self, built_in_filter) -> set:
        q = get_extra_filters_q(self.analysis, built_in_filter)
        return set(Variant.objects.filter(q).values_list("pk", flat=True))

    def test_classified_matches_any_origin(self):
        self.assertEqual(self._filter_variant_pks(BuiltInFilters.CLASSIFIED), set(self.variant_pks.values()))

    def test_classified_pathogenic_is_germline_only(self):
        self.assertEqual(self._filter_variant_pks(BuiltInFilters.CLASSIFIED_PATHOGENIC),
                         {self.variant_pks[self.germline_p.pk]})

    def test_classified_tier_1_2_is_somatic_only(self):
        self.assertEqual(self._filter_variant_pks(BuiltInFilters.CLASSIFIED_TIER_1_2),
                         {self.variant_pks[self.somatic_tier_1.pk],
                          self.variant_pks[self.somatic_tier_1_or_2.pk]})


class ClassificationsNodeInputTest(AnalysisSetupMixin, TestCase):
    def test_source_mode_takes_no_parent(self):
        node = ClassificationsNode(analysis=self.analysis, node_input=NodeMatchInput.MATCHING_VARIANTS)
        self.assertEqual((node.min_inputs, node.max_inputs), (0, 0))
        self.assertTrue(node.is_source)

    def test_filter_modes_take_one_parent(self):
        for node_input in (NodeMatchInput.PARENT_MATCHING, NodeMatchInput.PARENT_NOT_MATCHING):
            node = ClassificationsNode(analysis=self.analysis, node_input=node_input)
            self.assertEqual((node.min_inputs, node.max_inputs), (1, 1))
            self.assertFalse(node.is_source)

    def test_not_matching_negates_the_match(self):
        matching = ClassificationsNode.objects.create(analysis=self.analysis,
                                                      node_input=NodeMatchInput.PARENT_MATCHING)
        not_matching = ClassificationsNode.objects.create(analysis=self.analysis,
                                                          node_input=NodeMatchInput.PARENT_NOT_MATCHING)
        self.assertEqual(not_matching._get_node_q(), ~matching._get_node_q())


class ClassificationsNodeEditorTest(AnalysisSetupMixin, TestCase):
    """ The editor renders the allele origin toggle and one-line significance pill rows """

    def test_editor_renders(self):
        node = ClassificationsNode.objects.create(analysis=self.analysis)
        self.client.force_login(self.analysis.user)
        url = reverse('node_view', kwargs={"analysis_id": self.analysis.pk,
                                           "analysis_version": self.analysis.version,
                                           "node_id": node.pk,
                                           "node_version": node.version,
                                           "extra_filters": "default"})
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn("allele-origin-toggle", content)
        self.assertIn("germline-significances", content)
        self.assertIn("scs-tier_1", content)
        self.assertEqual(content.count('name="allele_origin"'), 3)
