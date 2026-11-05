"""
Coverage for the per-VAV backfill-completion flags that gate optimised query
branches in DamageNode. The flags live on VariantAnnotationVersion
(backfilled_spliceai_max_ds, backfilled_damage_counts, backfilled_ptc,
backfilled_sv_conservation) and are
flipped True by the matching `fix_*` / `backfill_*` management commands once a
partition's derived columns are populated.
"""
from django.test import TestCase, override_settings

from analysis.enums import PTCNMDFilter
from analysis.models.nodes.filters.conservation_node import ConservationNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import AnalysisSetupMixin, create_annotated_variants, node_keeps
from annotation.enums import NMDEscapeStatus
from annotation.models.models import VariantAnnotation


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


class _AnnotatedVariantsMixin(_BackfillFlagMixin):
    """ A variant per ANNOTATIONS entry, its VariantAnnotation row holding those columns """
    ANNOTATIONS: dict[str, dict] = {}

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        cls.variants = create_annotated_variants(vav, cls.ANNOTATIONS)

    def _filtered(self, node: DamageNode) -> set[str]:
        return node_keeps(node, self.variants)


def _spliceai(ag=0.0, al=0.0, dg=0.0, dl=0.0) -> dict:
    # dbscsnv is in the same OR, so give it a score below every splice_min used here
    return {"spliceai_pred_ds_ag": ag, "spliceai_pred_ds_al": al,
            "spliceai_pred_ds_dg": dg, "spliceai_pred_ds_dl": dl,
            "spliceai_max_ds": max(ag, al, dg, dl), "dbscsnv_ada_score": 0.0, "dbscsnv_rf_score": 0.0}


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeSpliceAIBackfillFlag(_AnnotatedVariantsMixin, TestCase):
    ANNOTATIONS = {
        "acceptor_gain": _spliceai(ag=0.6),
        "donor_loss_only": _spliceai(dl=0.6),
        "low": _spliceai(ag=0.2, al=0.2, dg=0.2, dl=0.2),
        "no_spliceai": {"dbscsnv_ada_score": 0.0, "dbscsnv_rf_score": 0.0},
    }

    def _unbackfilled(self):
        """ A partition annotated before spliceai_max_ds existed has it null on every row """
        self._set_vav_flag(backfilled_spliceai_max_ds=False)
        VariantAnnotation.objects.filter(variant__in=self.variants.values()).update(spliceai_max_ds=None)

    def test_backfilled_emits_single_max_ds_q(self):
        self._set_vav_flag(backfilled_spliceai_max_ds=True)
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        self.assertEqual({"acceptor_gain", "donor_loss_only"}, self._filtered(node))

    def test_unbackfilled_emits_per_field_ds_qs(self):
        self._unbackfilled()
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        self.assertEqual({"acceptor_gain", "donor_loss_only"}, self._filtered(node))

    def test_unbackfilled_required_no_allow_null_omits_per_field_isnull(self):
        self._unbackfilled()
        node = DamageNode(analysis=self.analysis, splice_min=0.5, splice_required=True, splice_allow_null=False)
        self.assertEqual({"acceptor_gain", "donor_loss_only"}, self._filtered(node))

    def test_unbackfilled_required_allow_null_true_includes_isnulls(self):
        self._unbackfilled()
        node = DamageNode(analysis=self.analysis, splice_min=0.5, splice_required=True, splice_allow_null=True)
        self.assertEqual({"acceptor_gain", "donor_loss_only", "no_spliceai"}, self._filtered(node))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDamageNodeDamageCountsBackfillFlag(_AnnotatedVariantsMixin, TestCase):
    ANNOTATIONS = {
        "three_pathogenic": {"predictions_num_pathogenic": 3},
        "one_pathogenic": {"predictions_num_pathogenic": 1},
    }

    def test_backfilled_emits_predictions_num_pathogenic_clause(self):
        self._set_vav_flag(backfilled_damage_counts=True)
        node = DamageNode(analysis=self.analysis, damage_predictions_min=2)
        self.assertEqual({"three_pathogenic"}, self._filtered(node))

    def test_unbackfilled_omits_predictions_num_pathogenic_clause(self):
        """ On a pre-backfill VAV the counts are 0 on every row, so the filter is skipped rather than
            excluding the whole partition """
        self._set_vav_flag(backfilled_damage_counts=False)
        VariantAnnotation.objects.filter(variant__in=self.variants.values()).update(predictions_num_pathogenic=0)
        node = DamageNode(analysis=self.analysis, damage_predictions_min=2)
        self.assertEqual({"three_pathogenic", "one_pathogenic"}, self._filtered(node))


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
            any("latest annotation version" in w for w in warnings),
            f"expected latest annotation version pointer, got: {warnings}",
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
class TestDamageNodePTCBackfillFlag(_AnnotatedVariantsMixin, TestCase):
    """ ptc_nmd (#579) is columns_version 5 only, so force the VAV up. """
    ANNOTATIONS = {
        "escaping": {"nmd_escape_status": NMDEscapeStatus.ESCAPING},
        "predicted_nmd": {"nmd_escape_status": NMDEscapeStatus.PREDICTED_NMD},
        "not_calculated": {"nmd_escape_status": None},
    }

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
        node = self._cv5_node(ptc_nmd=PTCNMDFilter.ESCAPING)
        self.assertEqual({"escaping"}, self._filtered(node))

    def test_unbackfilled_omits_nmd_escape_status_clause(self):
        """On pre-backfill VAVs nmd_escape_status is null on every row, so we skip the filter
        entirely rather than excluding the whole partition."""
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node(ptc_nmd=PTCNMDFilter.ESCAPING)
        self.assertEqual({"escaping", "predicted_nmd", "not_calculated"}, self._filtered(node))

    def test_ptc_warning_on_unbackfilled_vav(self):
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node(ptc_nmd=PTCNMDFilter.ESCAPING)
        warnings = node.get_warnings()
        self.assertTrue(
            any("PTC NMD" in w for w in warnings),
            f"expected PTC NMD warning, got: {warnings}",
        )
        self.assertTrue(
            any("latest annotation version" in w for w in warnings),
            f"expected latest annotation version pointer, got: {warnings}",
        )

    def test_no_ptc_warning_when_backfilled(self):
        self._set_vav_flag(backfilled_ptc=True)
        node = self._cv5_node(ptc_nmd=PTCNMDFilter.ESCAPING)
        self.assertEqual([], [w for w in node.get_warnings() if "PTC NMD" in w])

    def test_ptc_warning_silent_without_filter(self):
        self._set_vav_flag(backfilled_ptc=False)
        node = self._cv5_node()  # no ptc_nmd
        self.assertEqual([], [w for w in node.get_warnings() if "PTC NMD" in w])


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestConservationNodeSVBackfillFlag(_BackfillFlagMixin, TestCase):
    """ #1657: a node filtering on conservation scores warns while the SV values are known wrong. """

    def test_warns_when_not_backfilled(self):
        self._set_vav_flag(backfilled_sv_conservation=False)
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.5)
        warnings = node.get_warnings()
        self.assertTrue(
            any("latest annotation version" in w for w in warnings),
            f"expected latest annotation version pointer, got: {warnings}",
        )

    def test_no_warning_when_backfilled(self):
        self._set_vav_flag(backfilled_sv_conservation=True)
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.5)
        self.assertEqual([], [w for w in node.get_warnings() if "conservation scores" in w])

    def test_no_warning_when_not_filtering(self):
        self._set_vav_flag(backfilled_sv_conservation=False)
        node = ConservationNode(analysis=self.analysis)  # sliders at default - filters nothing
        self.assertEqual([], [w for w in node.get_warnings() if "conservation scores" in w])
