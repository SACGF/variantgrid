from django.test import TestCase

from analysis.enums import NodeMatchInput
from analysis.forms.forms_nodes import ClinVarNodeForm
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.tests.utils import AnalysisSetupMixin, node_keeps
from annotation.enums import ClinVarOncogenicity, ClinVarReviewStatus, Pathogenicity
from annotation.models.models_clinvar import ClinVar
from library.django_utils.django_partition import temporary_db_table
from snpdb.enums import AlleleOriginFilterDefault, SomaticClinicalSignificance
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class ClinVarNodeInputTest(AnalysisSetupMixin, TestCase):
    def test_source_mode_takes_no_parent(self):
        node = ClinVarNode(analysis=self.analysis, node_input=NodeMatchInput.MATCHING_VARIANTS)
        self.assertEqual((node.min_inputs, node.max_inputs), (0, 0))
        self.assertTrue(node.is_source)

    def test_filter_modes_take_one_parent(self):
        for node_input in (NodeMatchInput.PARENT_MATCHING, NodeMatchInput.PARENT_NOT_MATCHING):
            node = ClinVarNode(analysis=self.analysis, node_input=node_input)
            self.assertEqual((node.min_inputs, node.max_inputs), (1, 1))
            self.assertFalse(node.is_source)


class ClinVarNodeQTest(AnalysisSetupMixin, TestCase):
    """ Which of a set of named ClinVar records each node setting keeps - one record per value on each axis,
        each with the review status a stars test needs, plus a variant ClinVar has no record for """
    S = ClinVarReviewStatus.CRITERIA_PROVIDED_SINGLE_SUBMITTER  # 1 star
    E = ClinVarReviewStatus.REVIEWED_BY_EXPERT_PANEL  # 3 stars
    N = ClinVarReviewStatus.NO_ASSERTION_PROVIDED  # 0 stars
    RECORDS = {
        "pathogenic": {"highest_pathogenicity": Pathogenicity.PATHOGENIC, "review_status": E},
        "likely_pathogenic": {"highest_pathogenicity": Pathogenicity.LIKELY_PATHOGENIC, "review_status": S},
        "vus": {"highest_pathogenicity": Pathogenicity.UNCERTAIN, "review_status": S},
        "likely_benign": {"highest_pathogenicity": Pathogenicity.LIKELY_BENIGN, "review_status": S},
        "benign": {"highest_pathogenicity": Pathogenicity.BENIGN, "review_status": S},
        "tier_1": {"somatic_tier": SomaticClinicalSignificance.TIER_1, "somatic_review_status": E},
        # 1 star on its germline axis only
        "tier_1_or_2": {"somatic_tier": SomaticClinicalSignificance.TIER_1_OR_2, "somatic_review_status": N,
                        "review_status": S},
        "tier_3": {"somatic_tier": SomaticClinicalSignificance.TIER_3, "somatic_review_status": S},
        "oncogenic": {"highest_oncogenicity": ClinVarOncogenicity.ONCOGENIC, "oncogenic_review_status": E},
        "likely_oncogenic": {"highest_oncogenicity": ClinVarOncogenicity.LIKELY_ONCOGENIC,
                             "oncogenic_review_status": S},
        "oncogenicity_uncertain": {"highest_oncogenicity": ClinVarOncogenicity.UNCERTAIN,
                                   "oncogenic_review_status": S},
    }
    NO_GERMLINE_CALL = {"tier_1", "tier_1_or_2", "tier_3", "oncogenic", "likely_oncogenic", "oncogenicity_uncertain"}

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        clinvar_version = cls.analysis.annotation_version.clinvar_version
        cls.variants = {"no_record": slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)}
        with temporary_db_table(ClinVar, clinvar_version.get_partition_table()):
            for i, (name, fields) in enumerate(cls.RECORDS.items(), start=1):
                variant = slowly_create_test_variant("1", 1000 + i * 100, "A", "T", cls.grch37)
                ClinVar.objects.create(version=clinvar_version, variant=variant,
                                       clinvar_variation_id=i, clinvar_allele_id=i, **fields)
                cls.variants[name] = variant

    def _filtered(self, node: ClinVarNode) -> set[str]:
        return node_keeps(node, self.variants)

    def test_bare_node_is_has_a_clinvar_record(self):
        self.assertEqual(set(self.RECORDS), self._filtered(ClinVarNode(analysis=self.analysis)))

    def test_not_matching_negates_the_whole_match(self):
        """ 'Not in ClinVar' - everything the matching node drops, including variants ClinVar has no record for """
        not_matching = ClinVarNode(analysis=self.analysis, germline_benign=False,
                                   node_input=NodeMatchInput.PARENT_NOT_MATCHING)
        self.assertEqual({"benign", "no_record"}, self._filtered(not_matching))

    def test_deselecting_one_pill_filters_the_row(self):
        node = ClinVarNode(analysis=self.analysis, germline_benign=False, germline_likely_benign=False,
                           germline_other=False)
        self.assertEqual({"pathogenic", "likely_pathogenic", "vus"}, self._filtered(node))

    def test_germline_catch_all_pill_selects_no_call(self):
        node = ClinVarNode(analysis=self.analysis, germline_pathogenic=False, germline_likely_pathogenic=False,
                           germline_uncertain=False, germline_likely_benign=False, germline_benign=False)
        self.assertEqual(self.NO_GERMLINE_CALL, self._filtered(node))

    def test_oncogenicity_subset(self):
        node = ClinVarNode(analysis=self.analysis, oncogenicity_uncertain=False,
                           oncogenicity_likely_benign=False, oncogenicity_benign=False)
        self.assertEqual({"oncogenic", "likely_oncogenic"}, self._filtered(node))

    def test_somatic_tier_1_pulls_in_tier_1_or_2(self):
        """ Tier I/II records might be either, so they match when Tier I or II is selected """
        node = ClinVarNode(analysis=self.analysis, somatic_tier_2=False, somatic_tier_3=False,
                           somatic_tier_4=False)
        self.assertEqual({"tier_1", "tier_1_or_2"}, self._filtered(node))

    def test_somatic_tier_3_does_not_pull_in_tier_1_or_2(self):
        node = ClinVarNode(analysis=self.analysis, somatic_tier_1=False, somatic_tier_2=False,
                           somatic_tier_4=False)
        self.assertEqual({"tier_3"}, self._filtered(node))

    def test_narrowed_axes_match_a_record_either_one_keeps(self):
        node = ClinVarNode(analysis=self.analysis, germline_uncertain=False, germline_likely_benign=False,
                           germline_benign=False, germline_other=False, oncogenicity_uncertain=False,
                           oncogenicity_likely_benign=False, oncogenicity_benign=False)
        self.assertEqual({"pathogenic", "likely_pathogenic", "oncogenic", "likely_oncogenic"}, self._filtered(node))

    def test_axis_with_every_pill_on_adds_nothing_to_the_or(self):
        """ Germline P/LP with somatic untouched is still germline P/LP, not 'or any somatic record' """
        node = ClinVarNode(analysis=self.analysis, germline_uncertain=False, germline_likely_benign=False,
                           germline_benign=False, germline_other=False)
        self.assertEqual({"pathogenic", "likely_pathogenic"}, self._filtered(node))

    def test_allele_origin_gates_the_rows(self):
        """ Germline only - the somatic pills stop applying without being retouched """
        node = ClinVarNode(analysis=self.analysis, allele_origin=AlleleOriginFilterDefault.GERMLINE,
                           somatic_tier_1=False)
        self.assertEqual(set(self.RECORDS), self._filtered(node))

    def test_stars_alone_read_all_3_axes(self):
        """ Matches ClinVar.is_expert_panel_or_greater taking the max across the axes """
        node = ClinVarNode(analysis=self.analysis, stars_min=3)
        self.assertEqual({"pathogenic", "tier_1", "oncogenic"}, self._filtered(node))

    def test_stars_read_the_axis_being_filtered(self):
        node = ClinVarNode(analysis=self.analysis, stars_min=1, somatic_tier_2=False, somatic_tier_3=False,
                           somatic_tier_4=False)
        self.assertEqual({"tier_1"}, self._filtered(node))

    def test_variation_ids(self):
        vus_variation_id = list(self.RECORDS).index("vus") + 1
        node = ClinVarNode(analysis=self.analysis, variation_ids=[vus_variation_id])
        self.assertEqual({"vus"}, self._filtered(node))


class ClinVarNodeChipsTest(AnalysisSetupMixin, TestCase):
    """ The card reads what the node filters on, so two nodes with the same query look the same """
    GERMLINE_P_LP = {"germline_uncertain": False, "germline_likely_benign": False,
                     "germline_benign": False, "germline_other": False}

    def test_row_labels_lead_their_pills(self):
        node = ClinVarNode(analysis=self.analysis, somatic_tier_3=False, somatic_tier_4=False,
                           oncogenicity_uncertain=False, oncogenicity_likely_benign=False,
                           oncogenicity_benign=False, **self.GERMLINE_P_LP)
        chips = node.get_node_chips()
        self.assertEqual([c.text for c in chips], ["Germline", "P", "LP", "Somatic", "I", "II", "O", "LO"])
        self.assertEqual([c.text for c in chips if c.row_break], ["Germline", "Somatic"])

    def test_allele_origin_only_shows_up_where_it_filters(self):
        """ Germline-only with every germline pill on is the same query as showing all """
        show_all = ClinVarNode(analysis=self.analysis, **self.GERMLINE_P_LP)
        germline_only = ClinVarNode(analysis=self.analysis, allele_origin=AlleleOriginFilterDefault.GERMLINE,
                                    **self.GERMLINE_P_LP)
        self.assertEqual(show_all._get_node_q(), germline_only._get_node_q())
        self.assertEqual(show_all.get_node_chips(), germline_only.get_node_chips())

    def test_default_node_has_no_chips(self):
        self.assertEqual(ClinVarNode(analysis=self.analysis).get_node_chips(), [])


class ClinVarNodeFormVariationIdsTest(AnalysisSetupMixin, TestCase):
    def _clean_variation_ids(self, value):
        node = ClinVarNode.objects.create(analysis=self.analysis)
        form = ClinVarNodeForm(data={"allele_origin": AlleleOriginFilterDefault.SHOW_ALL,
                                     "node_input": NodeMatchInput.MATCHING_VARIANTS,
                                     "stars_min": 0, "variation_ids": value},
                               instance=node)
        form.is_valid()
        return form

    def test_separators(self):
        form = self._clean_variation_ids("12345, 67890\n111")
        self.assertEqual(form.cleaned_data["variation_ids"], [12345, 67890, 111])

    def test_blank(self):
        form = self._clean_variation_ids("")
        self.assertEqual(form.cleaned_data["variation_ids"], [])

    def test_non_numeric_rejected(self):
        form = self._clean_variation_ids("12345, banana")
        self.assertIn("variation_ids", form.errors)


class ClinVarNodeFormRememberedPillsTest(AnalysisSetupMixin, TestCase):
    """ The greyed-out pills for the origin the node isn't filtering on are disabled, so they don't post -
        the form has to keep the stored values rather than let them clear to False """

    def _save(self, node, allele_origin, **data):
        form = ClinVarNodeForm(data={"allele_origin": allele_origin,
                                     "node_input": NodeMatchInput.MATCHING_VARIANTS,
                                     "stars_min": 0, "variation_ids": "", **data},
                               instance=node)
        self.assertTrue(form.is_valid(), form.errors)
        return form.save()

    def test_germline_pills_survive_saving_in_somatic_mode(self):
        node = ClinVarNode.objects.create(analysis=self.analysis)
        self._save(node, AlleleOriginFilterDefault.SHOW_ALL,
                   germline_pathogenic="on", germline_likely_pathogenic="on",
                   somatic_tier_1="on", oncogenicity_oncogenic="on")
        node.refresh_from_db()
        self.assertFalse(node.germline_uncertain)

        # Somatic mode posts no germline pills at all - they should be left as they were
        self._save(node, AlleleOriginFilterDefault.SOMATIC, somatic_tier_2="on")
        node.refresh_from_db()
        self.assertTrue(node.germline_pathogenic)
        self.assertTrue(node.germline_likely_pathogenic)
        self.assertFalse(node.germline_uncertain)
        self.assertTrue(node.somatic_tier_2)
        self.assertFalse(node.somatic_tier_1)

    def test_somatic_pills_survive_saving_in_germline_mode(self):
        node = ClinVarNode.objects.create(analysis=self.analysis)
        self._save(node, AlleleOriginFilterDefault.SHOW_ALL,
                   germline_pathogenic="on", somatic_tier_1="on", oncogenicity_oncogenic="on")

        self._save(node, AlleleOriginFilterDefault.GERMLINE, germline_benign="on")
        node.refresh_from_db()
        self.assertTrue(node.somatic_tier_1)
        self.assertTrue(node.oncogenicity_oncogenic)
        self.assertFalse(node.oncogenicity_benign)
        self.assertTrue(node.germline_benign)
        self.assertFalse(node.germline_pathogenic)

    def test_show_all_saves_every_pill(self):
        node = ClinVarNode.objects.create(analysis=self.analysis)
        self._save(node, AlleleOriginFilterDefault.SHOW_ALL, germline_pathogenic="on", somatic_tier_1="on")
        node.refresh_from_db()
        self.assertTrue(node.germline_pathogenic)
        self.assertTrue(node.somatic_tier_1)
        self.assertFalse(node.oncogenicity_oncogenic)

    def test_matching_variants_label_names_the_source(self):
        form = ClinVarNodeForm(instance=ClinVarNode(analysis=self.analysis))
        choices = dict(form.fields["node_input"].choices)
        self.assertEqual(choices[NodeMatchInput.MATCHING_VARIANTS], "All records in ClinVar (no parent)")
