from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import create_annotated_variants, node_keeps
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import VariantAnnotationVersion
from annotation.pathogenicity_predictions import TOOLS, PathogenicityTool, RawScoreDirection
from snpdb.models.models_genome import GenomeBuild

RAW_SCORE_TOOLS = [tool for tool in TOOLS if tool.raw_field]
# (field, damaging value, other value)
PRED_FIELDS = [("alphamissense_pred", "P", "B"), ("clinpred_pred", "D", "T"),
               ("metarnn_pred", "D", "T"), ("primateai_pred", "D", "T")]


def _threshold(tool: PathogenicityTool) -> float:
    """ Three quarters of the way up the slider - inside every tool's range, and never 0 """
    return tool.raw_min + (tool.raw_max - tool.raw_min) * 0.75


def _raw_score(tool: PathogenicityTool, steps_to_damaging_side: int) -> float:
    """ A score the given number of slider steps from the threshold, towards the damaging side """
    offset = steps_to_damaging_side * tool.raw_step
    if tool.raw_direction == RawScoreDirection.LOWER:
        offset = -offset
    return _threshold(tool) + offset


class _ColumnsVersionMixin:
    """ columns_version is a derived property on AnalysisNode (reads through the VAV), so forcing
        the VAV is the only way to flip the node into a version's code path. """
    COLUMNS_VERSION = 4
    ROWS: dict[str, dict] = {}

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        vav = VariantAnnotationVersion.objects.filter(genome_build=cls.grch37).first()
        vav.columns_version = cls.COLUMNS_VERSION
        vav.save()
        user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        cls.variants = create_annotated_variants(vav, cls.ROWS)

    def _filtered(self, node: DamageNode, among: set[str]) -> set[str]:
        return node_keeps(node, self.variants, among)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV4Q(_ColumnsVersionMixin, TestCase):
    """ Verifies DamageNode raw-score and pred-field filtering at columns_version 4. """
    RAW_SCORE_ROWS = {"damaging", "at_threshold", "benign", "unscored"}
    PRED_ROWS = {"pred_damaging", "pred_other", "pred_null"}
    BAYESDEL_ROWS = {"bayesdel_negative", "bayesdel_zero", "bayesdel_positive"}
    ROWS = {
        "damaging": {tool.raw_field: _raw_score(tool, 1) for tool in RAW_SCORE_TOOLS},
        "at_threshold": {tool.raw_field: _raw_score(tool, 0) for tool in RAW_SCORE_TOOLS},
        "benign": {tool.raw_field: _raw_score(tool, -1) for tool in RAW_SCORE_TOOLS},
        "unscored": {},
        "pred_damaging": {field: damaging for field, damaging, _ in PRED_FIELDS},
        "pred_other": {field: other for field, _, other in PRED_FIELDS},
        "pred_null": {},
        # BayesDel's range spans 0, so a slider at 0 is a real threshold
        "bayesdel_negative": {"bayesdel_noaf_score": -0.5},
        "bayesdel_zero": {"bayesdel_noaf_score": 0.0},
        "bayesdel_positive": {"bayesdel_noaf_score": 0.5},
    }

    # Several tools are GRCh38 only, and a column the version doesn't populate isn't filtered on (#2787)
    @patch.object(VariantAnnotationVersion, "visible_columns", frozenset(t.raw_field for t in RAW_SCORE_TOOLS))
    def test_raw_score_threshold_keeps_damaging_side_per_tool(self):
        for tool in RAW_SCORE_TOOLS:
            with self.subTest(tool=tool.name):
                node = DamageNode(analysis=self.analysis, **{tool.node_threshold_field: _threshold(tool)})
                self.assertEqual({"damaging", "at_threshold"}, self._filtered(node, self.RAW_SCORE_ROWS))

    def test_raw_score_threshold_of_zero_filters(self):
        """ The editor's enable box clears a slider to None; an enabled slider at 0 is a threshold (#2657) """
        node = DamageNode(analysis=self.analysis, bayesdel_noaf_score_min=0.0)
        self.assertTrue(node.modifies_parents())
        self.assertEqual({"bayesdel_zero", "bayesdel_positive"}, self._filtered(node, self.BAYESDEL_ROWS))

    def test_raw_score_required_at_zero_keeps_only_scored(self):
        """ On a 0..1 scale, required >= 0 with allow-null off is "has a score": every scored variant passes """
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.0,
                          alphamissense_score_required=True, alphamissense_score_allow_null=False)
        self.assertEqual({"damaging", "at_threshold", "benign"}, self._filtered(node, self.RAW_SCORE_ROWS))

    def test_pred_field_keeps_matching_prediction_per_tool(self):
        for field, damaging, _ in PRED_FIELDS:
            with self.subTest(field=field):
                node = DamageNode(analysis=self.analysis, **{field: damaging})
                self.assertEqual({"pred_damaging"}, self._filtered(node, self.PRED_ROWS))

    def test_raw_score_required_with_allow_null_keeps_unscored(self):
        alphamissense = next(tool for tool in TOOLS if tool.raw_field == "alphamissense_score")
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=_threshold(alphamissense),
                          alphamissense_score_required=True, alphamissense_score_allow_null=True)
        self.assertEqual({"damaging", "at_threshold", "unscored"}, self._filtered(node, self.RAW_SCORE_ROWS))
        node.alphamissense_score_allow_null = False
        self.assertEqual({"damaging", "at_threshold"}, self._filtered(node, self.RAW_SCORE_ROWS))


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV3IgnoresV4Fields(_ColumnsVersionMixin, TestCase):
    """ A v3 (or earlier) DamageNode must ignore v4-only model fields, even if set. """
    COLUMNS_VERSION = 3
    ROWS = {"alphamissense_high": {"alphamissense_score": 0.9},
            "alphamissense_low": {"alphamissense_score": 0.05}}

    def test_v3_ignores_raw_score_threshold(self):
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.5)
        self.assertFalse(node.modifies_parents())
        self.assertEqual({"alphamissense_high", "alphamissense_low"}, self._filtered(node, set(self.ROWS)))


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeVersionWithoutPredictors(_ColumnsVersionMixin, TestCase):
    """ #2787: T2T has no dbNSFP, so every rankscore is null and every variant counts 0 predictions - a filter
        on them would drop the whole sample. Skipped with a warning instead, applied where the version has them. """
    COLUMNS_VERSION = 3
    ROWS = {"revel_high": {"revel_rankscore": 0.9, "predictions_num_pathogenic": 1},
            "revel_low": {"revel_rankscore": 0.1, "predictions_num_pathogenic": 0}}

    def _node(self):
        return DamageNode(analysis=self.analysis, revel_rankscore_min=0.5, damage_predictions_min=1)

    def test_applied_where_populated(self):
        with patch.object(VariantAnnotationVersion, "visible_columns", frozenset({"revel_rankscore"})):
            node = self._node()
            self.assertEqual({"revel_high"}, self._filtered(node, set(self.ROWS)))
            self.assertFalse([w for w in node.get_warnings() if "NOT being applied" in w])

    def test_skipped_with_warning_where_not_populated(self):
        with patch.object(VariantAnnotationVersion, "visible_columns", frozenset()):
            node = self._node()
            self.assertEqual(set(self.ROWS), self._filtered(node, set(self.ROWS)))
            self.assertEqual(2, len([w for w in node.get_warnings() if "NOT being applied" in w]))
            self.assertFalse(DamageNode(analysis=self.analysis).show_rankscore_sliders())
