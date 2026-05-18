from django.test import TestCase, override_settings

from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.tests.utils import AnalysisSetupMixin


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeMaxAfGate(AnalysisSetupMixin, TestCase):

    def _pop_node(self, **kwargs):
        defaults = dict(
            percent=1.0,
            gnomad_af=True,
            gnomad_popmax_af=False,
            af_1kg=True,
            af_uk10k=True,
            topmed_af=False,
        )
        defaults.update(kwargs)
        return PopulationNode.objects.create(analysis=self.analysis, **defaults)

    def test_filtering_by_population_emits_max_af_gate(self):
        node = self._pop_node(percent=1.0)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__max_af__isnull", q_str)
        self.assertIn("variantannotation__max_af__lte", q_str)
        self.assertIn("0.01", q_str)

    def test_filtering_disabled_does_not_emit_max_af_gate(self):
        node = self._pop_node(percent=PopulationNode.EVERYTHING)
        q = node._get_node_q()
        q_str = str(q) if q is not None else ""
        self.assertNotIn("variantannotation__max_af", q_str)

    def test_per_field_clauses_still_emitted(self):
        node = self._pop_node(percent=1.0)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__gnomad_af__lte", q_str)
        self.assertIn("variantannotation__af_1kg__lte", q_str)
        self.assertIn("variantannotation__af_uk10k__lte", q_str)

    def test_max_af_gate_uses_percent_over_100(self):
        node = self._pop_node(percent=5.0)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__max_af__lte", q_str)
        self.assertIn("0.05", q_str)

    def test_max_af_gate_skipped_on_unbackfilled_vav(self):
        """Pre-backfill VAVs (backfilled_max_af=False) skip the max_af gate clause
        but per-field clauses are unaffected. Profiling code can also flip the flag
        temporarily for A/B comparison."""
        vav = self.analysis.annotation_version.variant_annotation_version
        original = vav.backfilled_max_af
        vav.backfilled_max_af = False
        vav.save(update_fields=["backfilled_max_af"])
        try:
            node = self._pop_node(percent=1.0)
            q_str = str(node._get_node_q())
            self.assertNotIn("variantannotation__max_af__lte", q_str)
            self.assertNotIn("variantannotation__max_af__isnull", q_str)
            self.assertIn("variantannotation__gnomad_af__lte", q_str)
        finally:
            vav.backfilled_max_af = original
            vav.save(update_fields=["backfilled_max_af"])
