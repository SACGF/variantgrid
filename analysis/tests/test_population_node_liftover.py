"""
PopulationNode on gnomAD's own T2T liftover (#2788): variants where CHM13's reference allele differs from GRCh38's
have no gnomAD AF, so a gnomAD frequency filter keeps them as rare - the node warns, until the version uses the
regenerated file that adds them.
"""
from django.test import TestCase

from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.tests.utils import AnalysisSetupMixin

LIFTOVER_WARNING = "is lifted over from GRCh38"


class TestPopulationNodeGnomADLiftover(AnalysisSetupMixin, TestCase):
    def _set_gnomad(self, gnomad: str):
        vav = self.analysis.annotation_version.variant_annotation_version
        vav.gnomad = gnomad
        vav.save()

    def _warns(self, **kwargs) -> bool:
        node = PopulationNode(analysis=self.analysis, **kwargs)
        return any(LIFTOVER_WARNING in w for w in node.get_warnings())

    def test_warns_when_filtering_on_lifted_over_gnomad(self):
        self._set_gnomad("4.1.t2t_liftover")
        self.assertTrue(self._warns(percent=1))
        self.assertTrue(self._warns(gnomad_hom_alt_max=2))

    def test_no_warning_unless_filtering_on_gnomad(self):
        self._set_gnomad("4.1.t2t_liftover")
        self.assertFalse(self._warns())  # 100% - no filter
        self.assertFalse(self._warns(percent=1, gnomad_af=False, af_1kg=True))

    def test_no_warning_on_native_gnomad_or_with_ref_swaps(self):
        for gnomad in ["4.1", "4.1.t2t_liftover_ref_swaps"]:
            with self.subTest(gnomad=gnomad):
                self._set_gnomad(gnomad)
                self.assertFalse(self._warns(percent=1))
