from django.db.models import Q
from django.test import TestCase

from analysis.forms.forms_nodes import ClassificationsNodeForm
from analysis.models.enums import ClassificationsNodeInput, ClinVarRecordFilter
from analysis.models.nodes.sources.classifications_node import ClassificationsNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity, ClinVarReviewStatus
from classification.enums import SomaticClinicalSignificance


class ClassificationsNodeInputTest(AnalysisSetupMixin, TestCase):
    def test_source_mode_takes_no_parent(self):
        node = ClassificationsNode(analysis=self.analysis, node_input=ClassificationsNodeInput.MATCHING_VARIANTS)
        self.assertEqual((node.min_inputs, node.max_inputs), (0, 0))
        self.assertTrue(node.is_source)

    def test_filter_modes_take_one_parent(self):
        for node_input in (ClassificationsNodeInput.PARENT_MATCHING, ClassificationsNodeInput.PARENT_NOT_MATCHING):
            node = ClassificationsNode(analysis=self.analysis, node_input=node_input)
            self.assertEqual((node.min_inputs, node.max_inputs), (1, 1))
            self.assertFalse(node.is_source)

    def test_not_matching_negates_the_whole_match(self):
        """ 'Neither our labs nor ClinVar call this benign' - the negation wraps the OR, not one side """
        kwargs = {"analysis": self.analysis, "clinvar_benign": True}
        matching = ClassificationsNode.objects.create(node_input=ClassificationsNodeInput.PARENT_MATCHING, **kwargs)
        not_matching = ClassificationsNode.objects.create(node_input=ClassificationsNodeInput.PARENT_NOT_MATCHING,
                                                          **kwargs)
        self.assertEqual(not_matching._get_node_q(), ~matching._get_node_q())


class ClassificationsNodeClinVarQTest(AnalysisSetupMixin, TestCase):
    def test_no_clinvar_settings_contributes_nothing(self):
        node = ClassificationsNode.objects.create(analysis=self.analysis)
        self.assertFalse(node.has_clinvar_filters())
        self.assertEqual(node._clinvar_q(), Q())
        # An empty Q drops out of the OR, so the node still asks exactly what it asked before
        self.assertEqual(node._get_node_q(), node._classifications_q())

    def test_significance_include(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_pathogenic=True, clinvar_likely_pathogenic=True)
        self.assertEqual(node._clinvar_q(),
                         Q(clinvar__highest_pathogenicity__in=[ClinVarPathogenicity.LIKELY_PATHOGENIC,
                                                              ClinVarPathogenicity.PATHOGENIC]))

    def test_significance_exclude_is_a_negated_subquery(self):
        """ ~Q compiles to NOT EXISTS, so variants with no ClinVar record survive - unlike a `ne` lookup """
        node = ClassificationsNode(analysis=self.analysis, clinvar_benign=True, clinvar_likely_benign=True,
                                   clinvar_significance_exclude=True)
        self.assertEqual(node._clinvar_q(),
                         ~Q(clinvar__highest_pathogenicity__in=[ClinVarPathogenicity.BENIGN,
                                                                ClinVarPathogenicity.LIKELY_BENIGN]))

    def test_exclude_benign_paired_with_has_record(self):
        """ The faithful replacement for the prod rule 'highest_pathogenicity != 1 AND != 2' """
        node = ClassificationsNode(analysis=self.analysis, clinvar_benign=True, clinvar_likely_benign=True,
                                   clinvar_significance_exclude=True,
                                   clinvar_record=ClinVarRecordFilter.PRESENT)
        expected = ~Q(clinvar__highest_pathogenicity__in=[ClinVarPathogenicity.BENIGN,
                                                          ClinVarPathogenicity.LIKELY_BENIGN]) \
            & Q(clinvar__isnull=False)
        self.assertEqual(node._clinvar_q(), expected)

    def test_absent_record(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_record=ClinVarRecordFilter.ABSENT)
        self.assertEqual(node._clinvar_q(), Q(clinvar__isnull=True))

    def test_variation_ids(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_variation_ids=[12345])
        self.assertTrue(node.has_clinvar_filters())
        self.assertEqual(node._clinvar_q(), Q(clinvar__clinvar_variation_id__in=[12345]))


class ClassificationsNodeFormVariationIdsTest(AnalysisSetupMixin, TestCase):
    def _clean_variation_ids(self, value):
        node = ClassificationsNode.objects.create(analysis=self.analysis)
        form = ClassificationsNodeForm(data={"allele_origin": "A", "node_input": ClassificationsNodeInput.MATCHING_VARIANTS,
                                             "clinvar_record": ClinVarRecordFilter.ANY, "clinvar_stars_min": 0,
                                             "clinvar_variation_ids": value},
                                       instance=node)
        form.is_valid()
        return form

    def test_separators(self):
        form = self._clean_variation_ids("12345, 67890\n111")
        self.assertEqual(form.cleaned_data["clinvar_variation_ids"], [12345, 67890, 111])

    def test_blank(self):
        form = self._clean_variation_ids("")
        self.assertEqual(form.cleaned_data["clinvar_variation_ids"], [])

    def test_non_numeric_rejected(self):
        form = self._clean_variation_ids("12345, banana")
        self.assertIn("clinvar_variation_ids", form.errors)


class ClassificationsNodeClinVarSomaticQTest(AnalysisSetupMixin, TestCase):
    def test_somatic_tier(self):
        """ Tier I/II records might be either, so they match when Tier I or II is selected """
        node = ClassificationsNode(analysis=self.analysis, clinvar_tier_1=True)
        self.assertTrue(node.has_clinvar_filters())
        self.assertEqual(node._clinvar_q(),
                         Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_1,
                                                      SomaticClinicalSignificance.TIER_1_OR_2]))

    def test_somatic_tier_3_does_not_pull_in_tier_1_or_2(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_tier_3=True)
        self.assertEqual(node._clinvar_q(),
                         Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_3]))

    def test_oncogenicity(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_likely_oncogenic=True, clinvar_oncogenic=True)
        self.assertEqual(node._clinvar_q(),
                         Q(clinvar__highest_oncogenicity__in=[ClinVarOncogenicity.LIKELY_ONCOGENIC,
                                                             ClinVarOncogenicity.ONCOGENIC]))

    def test_stars_read_the_axis_being_filtered(self):
        node = ClassificationsNode(analysis=self.analysis, clinvar_tier_1=True, clinvar_stars_min=1)
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(1)
        expected = Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_1,
                                                SomaticClinicalSignificance.TIER_1_OR_2]) \
            & Q(clinvar__somatic_review_status__in=review_statuses)
        self.assertEqual(node._clinvar_q(), expected)

    def test_stars_alone_read_all_3_axes(self):
        """ Matches ClinVar.is_expert_panel_or_greater taking the max across the axes """
        node = ClassificationsNode(analysis=self.analysis, clinvar_stars_min=3)
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(3)
        expected = Q(clinvar__review_status__in=review_statuses) \
            | Q(clinvar__somatic_review_status__in=review_statuses) \
            | Q(clinvar__oncogenic_review_status__in=review_statuses)
        self.assertEqual(node._clinvar_q(), expected)
