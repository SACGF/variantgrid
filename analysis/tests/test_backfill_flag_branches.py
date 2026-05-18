"""
Coverage for the per-VAV backfill-completion flags that gate optimised query
branches in DamageNode and PopulationNode. The flags live on
VariantAnnotationVersion (backfilled_spliceai_max_ds, backfilled_max_af,
backfilled_damage_counts) and are flipped True by the matching `fix_*`
management commands once a partition's derived columns are populated.
"""
from django.test import TestCase, override_settings

from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.tests.utils import AnalysisSetupMixin


class _BackfillFlagMixin(AnalysisSetupMixin):
    """Helper to flip a VAV flag for the duration of a single test and restore it."""

    def _set_vav_flag(self, **flags):
        """Bulk-update flag fields on this analysis's VAV. Returns a restore callable."""
        vav = self.analysis.annotation_version.variant_annotation_version
        original = {name: getattr(vav, name) for name in flags}
        for name, value in flags.items():
            setattr(vav, name, value)
        vav.save(update_fields=list(flags.keys()))

        def restore():
            for name, value in original.items():
                setattr(vav, name, value)
            vav.save(update_fields=list(original.keys()))
        self.addCleanup(restore)
        return vav


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeSpliceAIBackfillFlag(_BackfillFlagMixin, TestCase):

    def test_backfilled_emits_single_max_ds_q(self):
        """backfilled_spliceai_max_ds=True → use the optimised spliceai_max_ds column."""
        self._set_vav_flag(backfilled_spliceai_max_ds=True)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__spliceai_max_ds__gte", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_ag", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_al", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_dg", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_dl", q_str)

    def test_unbackfilled_emits_per_field_ds_qs(self):
        """backfilled_spliceai_max_ds=False → fall back to per-DS-field loop."""
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        q_str = str(node._get_node_q())
        self.assertNotIn("variantannotation__spliceai_max_ds__gte", q_str)
        self.assertIn("variantannotation__spliceai_pred_ds_ag__gte", q_str)
        self.assertIn("variantannotation__spliceai_pred_ds_al__gte", q_str)
        self.assertIn("variantannotation__spliceai_pred_ds_dg__gte", q_str)
        self.assertIn("variantannotation__spliceai_pred_ds_dl__gte", q_str)

    def test_unbackfilled_required_allow_null_includes_per_field_isnull(self):
        """splice_required + splice_allow_null on pre-backfill VAV emits per-field isnull
        alternatives so rows with one qualifying DS field still match."""
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        node = DamageNode(
            analysis=self.analysis,
            splice_min=0.5,
            splice_required=True,
            splice_allow_null=False,
        )
        q_str = str(node._get_node_q())
        # allow_null=False → no per-field isnull alternative
        self.assertNotIn("spliceai_pred_ds_ag__isnull", q_str)
        # but the per-field gte clauses still emitted
        self.assertIn("variantannotation__spliceai_pred_ds_ag__gte", q_str)
        self.assertIn("variantannotation__spliceai_pred_ds_dl__gte", q_str)

    def test_unbackfilled_required_allow_null_true_includes_isnulls(self):
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        node = DamageNode(
            analysis=self.analysis,
            splice_min=0.5,
            splice_required=True,
            splice_allow_null=True,
        )
        q_str = str(node._get_node_q())
        for ds in ("spliceai_pred_ds_ag", "spliceai_pred_ds_al",
                   "spliceai_pred_ds_dg", "spliceai_pred_ds_dl"):
            self.assertIn(f"variantannotation__{ds}__isnull", q_str)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeDamageCountsBackfillFlag(_BackfillFlagMixin, TestCase):

    def test_backfilled_emits_predictions_num_pathogenic_clause(self):
        self._set_vav_flag(backfilled_damage_counts=True)
        node = DamageNode(analysis=self.analysis, damage_predictions_min=2)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__predictions_num_pathogenic__gte", q_str)

    def test_unbackfilled_omits_predictions_num_pathogenic_clause(self):
        """On pre-backfill VAV, predictions_num_pathogenic / predictions_num_benign
        default to 0 on every row, so we skip the filter entirely rather than
        excluding the whole partition."""
        self._set_vav_flag(backfilled_damage_counts=False)
        node = DamageNode(analysis=self.analysis, damage_predictions_min=2)
        q = node._get_node_q()
        q_str = str(q) if q is not None else ""
        self.assertNotIn("predictions_num_pathogenic", q_str)
        self.assertNotIn("predictions_num_benign", q_str)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeBackfillFlag(_BackfillFlagMixin, TestCase):

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

    def test_backfilled_emits_max_af_gate(self):
        self._set_vav_flag(backfilled_max_af=True)
        node = self._pop_node()
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__max_af__lte", q_str)
        self.assertIn("variantannotation__max_af__isnull", q_str)

    def test_unbackfilled_omits_max_af_gate(self):
        self._set_vav_flag(backfilled_max_af=False)
        node = self._pop_node()
        q_str = str(node._get_node_q())
        self.assertNotIn("variantannotation__max_af__lte", q_str)
        self.assertNotIn("variantannotation__max_af__isnull", q_str)
        # per-field clauses unaffected
        self.assertIn("variantannotation__gnomad_af__lte", q_str)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeBackfillFlagWarnings(_BackfillFlagMixin, TestCase):

    def test_no_warnings_when_backfilled(self):
        self._set_vav_flag(
            backfilled_spliceai_max_ds=True,
            backfilled_damage_counts=True,
        )
        node = DamageNode(
            analysis=self.analysis,
            splice_min=0.5,
            damage_predictions_min=2,
        )
        warnings = node.get_warnings()
        self.assertEqual([], [w for w in warnings if "spliceai" in w.lower() or "damage_predictions" in w.lower()])

    def test_damage_predictions_warning_on_unbackfilled_vav(self):
        self._set_vav_flag(backfilled_damage_counts=False)
        node = DamageNode(analysis=self.analysis, damage_predictions_min=2)
        warnings = node.get_warnings()
        self.assertTrue(
            any("damage_predictions_min" in w for w in warnings),
            f"expected damage_predictions_min warning, got: {warnings}",
        )
        self.assertTrue(
            any("fix_columns_version" in w for w in warnings),
            f"expected fix_columns_version pointer, got: {warnings}",
        )

    def test_damage_predictions_warning_silent_without_filter(self):
        """The warning only fires when the user has actually configured
        damage_predictions_min, so unrelated DamageNodes stay quiet."""
        self._set_vav_flag(backfilled_damage_counts=False)
        node = DamageNode(analysis=self.analysis)  # no damage_predictions_min
        warnings = node.get_warnings()
        self.assertFalse(
            any("damage_predictions_min" in w for w in warnings),
            f"warning fired without filter configured, got: {warnings}",
        )

    def test_spliceai_warning_on_unbackfilled_vav(self):
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        warnings = node.get_warnings()
        self.assertTrue(
            any("SpliceAI" in w for w in warnings),
            f"expected SpliceAI warning, got: {warnings}",
        )
        self.assertTrue(
            any("fix_historical_spliceai_max_ds" in w for w in warnings),
            f"expected fix_historical_spliceai_max_ds pointer, got: {warnings}",
        )

    def test_spliceai_warning_silent_without_filter(self):
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        node = DamageNode(analysis=self.analysis)  # no splice_min
        warnings = node.get_warnings()
        self.assertFalse(
            any("SpliceAI" in w for w in warnings),
            f"warning fired without filter configured, got: {warnings}",
        )
