""" The Effect node in words (#2836): the sentence its editor, summary and hover card show, built from the filters'
    groups - every "All of" filter and at least one "Any of" - which is the shape DamageNode._get_node_q runs """
from django.test import TestCase

from analysis.models import DamageNode
from analysis.models.nodes.filters.damage_node_filters import SENTENCE
from analysis.tests.utils import AnalysisSetupMixin
from annotation.enums import Pathogenicity, PathogenicityImpact
from snpdb.variant_filters import VariantType


class DamageNodeDescriptionTest(AnalysisSetupMixin, TestCase):

    def _describe(self, **kwargs) -> str:
        return DamageNode(analysis=self.analysis, **kwargs).get_node_description()

    def test_empty(self):
        self.assertEqual(SENTENCE["none"], self._describe())

    def test_any_of_only(self):
        self.assertEqual("Keeps variants with impact ≥ HIGH.", self._describe(impact_min=PathogenicityImpact.HIGH))
        self.assertEqual("Keeps variants with at least one of: impact ≥ MODERATE · REVEL ≥ 0.644.",
                         self._describe(impact_min=PathogenicityImpact.MODERATE, revel_score_min=0.644))

    def test_all_of_only(self):
        """ Variant type always narrows, so it is in All of with no _required of its own """
        self.assertEqual("Keeps variants with all of: impact ≥ HIGH · variant type is SNV.",
                         self._describe(impact_min=PathogenicityImpact.HIGH, impact_required=True,
                                        variant_types=[VariantType.SNV]))

    def test_pvs1_shape(self):
        description = self._describe(impact_min=PathogenicityImpact.HIGH,
                                     annotsv_acmg_class_min=Pathogenicity.LIKELY_PATHOGENIC,
                                     gene_loeuf_max=0.6, gene_loeuf_required=True)
        self.assertEqual("Keeps variants with at least one of: impact ≥ HIGH · AnnotSV ACMG class ≥ Likely "
                         "pathogenic, and gnomAD LOEUF ≤ 0.6 (or gene has no score).", description)

    def test_not_scored_only_in_all_of(self):
        """ _allow_null is True by default but only has an effect on a required filter """
        self.assertEqual("Keeps variants with REVEL ≥ 0.5.", self._describe(revel_score_min=0.5))
        self.assertEqual("Keeps variants with REVEL ≥ 0.5 (or not scored).",
                         self._describe(revel_score_min=0.5, revel_score_required=True))
        self.assertEqual("Keeps variants with REVEL ≥ 0.5.",
                         self._describe(revel_score_min=0.5, revel_score_required=True, revel_score_allow_null=False))

    def test_splicing_reads_both_thresholds(self):
        self.assertEqual("Keeps variants with splicing ≥ 0.5 or AlphaGenome splicing ≥ 1.",
                         self._describe(splice_min=0.5, alphagenome_splicing_min=1.0))

    def test_summary_is_the_description_in_html(self):
        node = DamageNode(analysis=self.analysis, impact_min=PathogenicityImpact.HIGH, gene_mis_z_min=3.09,
                          gene_mis_z_required=True, gene_mis_z_allow_null=False, in_cosmic=True)
        self.assertEqual("Keeps variants with <b>at least one of</b>: impact ≥ HIGH · in COSMIC, <b>and</b> "
                         "gnomAD missense Z ≥ 3.09.", node._get_method_summary())
