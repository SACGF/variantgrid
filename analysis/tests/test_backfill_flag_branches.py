"""
Coverage for the per-VAV backfill-completion flags that gate optimised query
branches in DamageNode. The flags live on VariantAnnotationVersion
(backfilled_spliceai_max_ds, backfilled_damage_counts, backfilled_ptc,
backfilled_sv_conservation) and are
flipped True by the matching `fix_*` / `backfill_*` management commands once a
partition's derived columns are populated.
"""
from django.test import TestCase, override_settings

from analysis.models.nodes.filters.conservation_node import ConservationNode
from analysis.models.nodes.filters.damage_node import DamageNode
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

    def test_unbackfilled_required_no_allow_null_omits_per_field_isnull(self):
        """splice_required without splice_allow_null on a pre-backfill VAV keeps the per-field
        gte clauses but emits no isnull alternatives."""
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


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodePTCBackfillFlag(_BackfillFlagMixin, TestCase):
    """ ptc_nmd_escaping (#579) is columns_version 5 only, so force the VAV up. """

    def _cv5_node(self, **kwargs) -> DamageNode:
        vav = self.analysis.annotation_version.variant_annotation_version
        original_columns_version = vav.columns_version
        vav.columns_version = 5
        vav.save(update_fields=["columns_version"])

        def restore():
            vav.columns_version = original_columns_version
            vav.save(update_fields=["columns_version"])
        self.addCleanup(restore)
        return DamageNode(analysis=self.analysis, **kwargs)

    def test_backfilled_emits_nmd_escape_status_clause(self):
        self._set_vav_flag(backfilled_ptc=True)
        node = self._cv5_node(ptc_nmd_escaping=True)
        self.assertIn("variantannotation__nmd_escape_status", str(node._get_node_q()))

    def test_unbackfilled_omits_nmd_escape_status_clause(self):
        """On pre-backfill VAVs nmd_escape_status is null on every row, so we skip the filter
        entirely rather than excluding the whole partition."""
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node(ptc_nmd_escaping=True)
        q = node._get_node_q()
        self.assertNotIn("nmd_escape_status", str(q) if q is not None else "")

    def test_ptc_warning_on_unbackfilled_vav(self):
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node(ptc_nmd_escaping=True)
        warnings = node.get_warnings()
        self.assertTrue(
            any("ptc_nmd_escaping" in w for w in warnings),
            f"expected ptc_nmd_escaping warning, got: {warnings}",
        )
        self.assertTrue(
            any("backfill_ptc_annotation" in w for w in warnings),
            f"expected backfill_ptc_annotation pointer, got: {warnings}",
        )

    def test_no_ptc_warning_when_backfilled(self):
        self._set_vav_flag(backfilled_ptc=True)
        node = self._cv5_node(ptc_nmd_escaping=True)
        self.assertEqual([], [w for w in node.get_warnings() if "ptc_nmd_escaping" in w])

    def test_ptc_warning_silent_without_filter(self):
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node()  # no ptc_nmd_escaping
        self.assertEqual([], [w for w in node.get_warnings() if "ptc_nmd_escaping" in w])


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestConservationNodeSVBackfillFlag(_BackfillFlagMixin, TestCase):
    """ #1657: a node filtering on conservation scores warns while the SV values are known wrong. """

    def test_warns_when_not_backfilled(self):
        self._set_vav_flag(backfilled_sv_conservation=False)
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.5)
        warnings = node.get_warnings()
        self.assertTrue(
            any("backfill_sv_conservation" in w for w in warnings),
            f"expected backfill_sv_conservation pointer, got: {warnings}",
        )

    def test_no_warning_when_backfilled(self):
        self._set_vav_flag(backfilled_sv_conservation=True)
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.5)
        self.assertEqual([], [w for w in node.get_warnings() if "backfill_sv_conservation" in w])

    def test_no_warning_when_not_filtering(self):
        self._set_vav_flag(backfilled_sv_conservation=False)
        node = ConservationNode(analysis=self.analysis)  # sliders at default - filters nothing
        self.assertEqual([], [w for w in node.get_warnings() if "backfill_sv_conservation" in w])
