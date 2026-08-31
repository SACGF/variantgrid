from django.db.models import Q
from django.test import TestCase

from analysis.forms.forms_nodes import ClinVarNodeForm
from analysis.models.enums import NodeMatchInput
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity, ClinVarReviewStatus
from classification.enums import SomaticClinicalSignificance
from snpdb.models.models_enums import AlleleOriginFilterDefault


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

    def test_bare_node_is_has_a_clinvar_record(self):
        node = ClinVarNode(analysis=self.analysis)
        self.assertEqual(node._get_node_q(), Q(clinvar__isnull=False))

    def test_not_matching_negates_the_whole_match(self):
        """ 'Not in ClinVar' - the negation wraps everything the node matched """
        matching = ClinVarNode(analysis=self.analysis, germline_benign=False)
        not_matching = ClinVarNode(analysis=self.analysis, germline_benign=False,
                                   node_input=NodeMatchInput.PARENT_NOT_MATCHING)
        self.assertEqual(not_matching._get_node_q(), ~matching._get_node_q())


class ClinVarNodeQTest(AnalysisSetupMixin, TestCase):
    HAS_RECORD = Q(clinvar__isnull=False)

    def test_deselecting_one_pill_filters_the_row(self):
        node = ClinVarNode(analysis=self.analysis, germline_benign=False, germline_likely_benign=False,
                           germline_other=False)
        expected = self.HAS_RECORD & Q(clinvar__highest_pathogenicity__in=[ClinVarPathogenicity.PATHOGENIC,
                                                                          ClinVarPathogenicity.LIKELY_PATHOGENIC,
                                                                          ClinVarPathogenicity.UNCERTAIN])
        self.assertEqual(node._get_node_q(), expected)

    def test_germline_catch_all_pill_selects_no_call(self):
        node = ClinVarNode(analysis=self.analysis, germline_pathogenic=False, germline_likely_pathogenic=False,
                           germline_uncertain=False, germline_likely_benign=False, germline_benign=False)
        self.assertEqual(node._get_node_q(), self.HAS_RECORD & Q(clinvar__highest_pathogenicity__in=[0]))

    def test_oncogenicity_catch_all_pill_selects_no_call(self):
        node = ClinVarNode(analysis=self.analysis, oncogenicity_oncogenic=False,
                           oncogenicity_likely_oncogenic=False, oncogenicity_uncertain=False,
                           oncogenicity_likely_benign=False, oncogenicity_benign=False)
        self.assertEqual(node._get_node_q(), self.HAS_RECORD & Q(clinvar__highest_oncogenicity__in=[0]))

    def test_somatic_catch_all_pill_selects_null_tier(self):
        node = ClinVarNode(analysis=self.analysis, somatic_tier_1=False, somatic_tier_2=False,
                           somatic_tier_3=False, somatic_tier_4=False)
        self.assertEqual(node._get_node_q(), self.HAS_RECORD & Q(clinvar__somatic_tier__isnull=True))

    def test_oncogenicity_subset(self):
        node = ClinVarNode(analysis=self.analysis, oncogenicity_uncertain=False,
                           oncogenicity_likely_benign=False, oncogenicity_benign=False,
                           oncogenicity_none=False)
        expected = self.HAS_RECORD & Q(clinvar__highest_oncogenicity__in=[ClinVarOncogenicity.ONCOGENIC,
                                                                          ClinVarOncogenicity.LIKELY_ONCOGENIC])
        self.assertEqual(node._get_node_q(), expected)

    def test_somatic_tier_1_pulls_in_tier_1_or_2(self):
        """ Tier I/II records might be either, so they match when Tier I or II is selected """
        node = ClinVarNode(analysis=self.analysis, somatic_tier_2=False, somatic_tier_3=False,
                           somatic_tier_4=False, somatic_tier_none=False)
        expected = self.HAS_RECORD & Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_1,
                                                                  SomaticClinicalSignificance.TIER_1_OR_2])
        self.assertEqual(node._get_node_q(), expected)

    def test_somatic_tier_3_does_not_pull_in_tier_1_or_2(self):
        node = ClinVarNode(analysis=self.analysis, somatic_tier_1=False, somatic_tier_2=False,
                           somatic_tier_4=False, somatic_tier_none=False)
        expected = self.HAS_RECORD & Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_3])
        self.assertEqual(node._get_node_q(), expected)

    def test_allele_origin_gates_the_rows(self):
        """ Germline only - the somatic pills stop applying without being retouched """
        node = ClinVarNode(analysis=self.analysis, allele_origin=AlleleOriginFilterDefault.GERMLINE,
                           somatic_tier_1=False)
        self.assertEqual(node._get_node_q(), self.HAS_RECORD)

    def test_stars_alone_read_all_3_axes(self):
        """ Matches ClinVar.is_expert_panel_or_greater taking the max across the axes """
        node = ClinVarNode(analysis=self.analysis, stars_min=3)
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(3)
        expected = self.HAS_RECORD & (Q(clinvar__review_status__in=review_statuses)
                                      | Q(clinvar__somatic_review_status__in=review_statuses)
                                      | Q(clinvar__oncogenic_review_status__in=review_statuses))
        self.assertEqual(node._get_node_q(), expected)

    def test_stars_read_the_axis_being_filtered(self):
        node = ClinVarNode(analysis=self.analysis, stars_min=1, somatic_tier_2=False, somatic_tier_3=False,
                           somatic_tier_4=False, somatic_tier_none=False)
        review_statuses = ClinVarReviewStatus.statuses_gte_stars(1)
        expected = self.HAS_RECORD \
            & Q(clinvar__somatic_tier__in=[SomaticClinicalSignificance.TIER_1,
                                           SomaticClinicalSignificance.TIER_1_OR_2]) \
            & Q(clinvar__somatic_review_status__in=review_statuses)
        self.assertEqual(node._get_node_q(), expected)

    def test_variation_ids(self):
        node = ClinVarNode(analysis=self.analysis, variation_ids=[12345])
        self.assertEqual(node._get_node_q(), self.HAS_RECORD & Q(clinvar__clinvar_variation_id__in=[12345]))


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
