"""
SACGF/variantgrid_sapath#410 — DamageNode surfaces a warning when its
splice filter runs against a VariantAnnotationVersion whose SpliceAI pin
is the raw precomputed file. Illumina's SpliceAI README recommends the
masked file for variant interpretation because the raw file also reports
strengthening of annotated splice sites and weakening of unannotated
ones, which are typically much less pathogenic.
"""
from django.test import TestCase, override_settings

from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import AnalysisSetupMixin


class _SpliceAIPinMixin(AnalysisSetupMixin):
    def _set_spliceai_pin(self, value: str):
        vav = self.analysis.annotation_version.variant_annotation_version
        original = vav.spliceai
        vav.spliceai = value
        vav.save(update_fields=["spliceai"])

        def restore():
            vav.spliceai = original
            vav.save(update_fields=["spliceai"])
        self.addCleanup(restore)
        return vav


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeRawSpliceAIWarning(_SpliceAIPinMixin, TestCase):

    def test_raw_spliceai_warning_when_splice_min_set(self):
        vav = self._set_spliceai_pin("raw 1.3")
        self.assertTrue(vav.uses_raw_spliceai)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        warnings = node.get_warnings()
        self.assertTrue(
            any("raw" in w and "SpliceAI" in w for w in warnings),
            f"expected raw SpliceAI warning, got: {warnings}",
        )

    def test_raw_spliceai_warning_silent_without_filter(self):
        self._set_spliceai_pin("raw 1.3")
        node = DamageNode(analysis=self.analysis)  # no splice_min
        warnings = node.get_warnings()
        self.assertFalse(
            any("raw" in w and "SpliceAI" in w for w in warnings),
            f"warning fired without filter configured, got: {warnings}",
        )

    def test_masked_spliceai_no_warning(self):
        vav = self._set_spliceai_pin("masked 1.3")
        self.assertFalse(vav.uses_raw_spliceai)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        warnings = node.get_warnings()
        self.assertFalse(
            any("raw" in w and "SpliceAI" in w for w in warnings),
            f"warning fired for masked VAV, got: {warnings}",
        )
