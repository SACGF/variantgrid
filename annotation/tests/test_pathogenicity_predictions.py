import json

from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from annotation.pathogenicity_predictions import (
    TOOLS,
    pred_pathogenic_funcs,
    raw_score_pathogenic_funcs,
)
from annotation.templatetags.pathogenicity_tags import pathogenicity_thresholds
from snpdb.models.models_genome import GenomeBuild


class TestPathogenicityPredictionsTable(TestCase):
    def test_raw_score_funcs_match_tools_with_thresholds(self):
        funcs = raw_score_pathogenic_funcs()
        for tool in TOOLS:
            if tool.raw_field and tool.raw_pathogenic_threshold is not None:
                self.assertIn(tool.raw_field, funcs)
                # Above threshold → pathogenic, below → not.
                self.assertTrue(funcs[tool.raw_field](tool.raw_pathogenic_threshold))
                self.assertFalse(funcs[tool.raw_field](tool.raw_pathogenic_threshold - 0.001))

    def test_varity_er_excluded_from_raw_score_funcs(self):
        # VARITY_ER has no calibrated threshold and must not contribute to the count.
        funcs = raw_score_pathogenic_funcs()
        self.assertNotIn("varity_er_score", funcs)

    def test_pred_funcs_count_alphamissense_lp_and_p(self):
        funcs = pred_pathogenic_funcs()
        self.assertIn("alphamissense_pred", funcs)
        self.assertTrue(funcs["alphamissense_pred"]("P"))
        self.assertTrue(funcs["alphamissense_pred"]("p"))  # likely_pathogenic
        self.assertFalse(funcs["alphamissense_pred"]("B"))
        self.assertFalse(funcs["alphamissense_pred"]("a"))  # ambiguous

    def test_pred_funcs_categorical_damaging(self):
        funcs = pred_pathogenic_funcs()
        for f in ("clinpred_pred", "metarnn_pred", "primateai_pred"):
            self.assertTrue(funcs[f]("D"))
            self.assertFalse(funcs[f]("T"))


class TestVariantAnnotationVersionRawScoreFuncs(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.vav = VariantAnnotationVersion.objects.filter(genome_build=cls.grch37).first()

    def test_returns_empty_pre_v4(self):
        for cv in (1, 2, 3):
            self.vav.columns_version = cv
            self.assertEqual(self.vav.get_raw_score_pathogenic_prediction_funcs(), {})

    def test_returns_raw_and_pred_funcs_at_v4(self):
        self.vav.columns_version = 4
        funcs = self.vav.get_raw_score_pathogenic_prediction_funcs()
        # Calibrated raw fields
        for raw_field in ("alphamissense_score", "bayesdel_noaf_score", "cadd_phred",
                          "clinpred_score", "metarnn_score", "mpc_score", "mutpred2_score",
                          "primateai_score", "revel_score", "varity_r_score", "vest4_score"):
            self.assertIn(raw_field, funcs)
        # Pred fields
        for pred_field in ("alphamissense_pred", "clinpred_pred",
                           "metarnn_pred", "primateai_pred"):
            self.assertIn(pred_field, funcs)

    def test_rankscore_funcs_unchanged_at_v2_v3(self):
        self.vav.columns_version = 2
        v2 = self.vav.get_rankscore_pathogenic_prediction_funcs()
        self.assertIn("revel_rankscore", v2)
        self.assertNotIn("alphamissense_rankscore", v2)
        self.vav.columns_version = 3
        v3 = self.vav.get_rankscore_pathogenic_prediction_funcs()
        self.assertIn("alphamissense_rankscore", v3)


class TestPathogenicityThresholdsTag(TestCase):
    def test_includes_pejaver_calibrated_revel_band(self):
        bands = json.loads(pathogenicity_thresholds())
        self.assertAlmostEqual(bands["revel_score"][0], 0.290)
        self.assertAlmostEqual(bands["revel_score"][1], 0.644)

    def test_includes_pejaver_calibrated_cadd_phred_band(self):
        bands = json.loads(pathogenicity_thresholds())
        self.assertEqual(bands["cadd_phred"], [22.7, 25.3])

    def test_excludes_uncalibrated_tools(self):
        bands = json.loads(pathogenicity_thresholds())
        self.assertNotIn("varity_er_score", bands)
