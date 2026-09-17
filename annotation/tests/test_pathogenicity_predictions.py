import json

from django.conf import settings
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
    """ The tag is variant_details.html CODE_THRESHOLDS - every coloured number on the page """

    def test_excludes_uncalibrated_tools(self):
        # A tool with no calibrated band must not get one
        bands = json.loads(pathogenicity_thresholds())
        self.assertNotIn("varity_er_score", bands)
        self.assertNotIn("phylop_46_way_mammalian", bands)  # COLOUR_BANDS entry, no cutoff chosen

    def test_rankscores_banded_off_settings(self):
        bands = json.loads(pathogenicity_thresholds())
        for rankscore_field in (t.rankscore_field for t in TOOLS if t.rankscore_field):
            self.assertEqual(bands[rankscore_field],
                             {"pathogenic": settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE,
                              "benign": settings.ANNOTATION_MAX_BENIGN_RANKSCORE})

    def test_rankscore_bands_cover_damage_count_columns(self):
        # The banded rankscores and the ones counted into predictions_num_pathogenic are the same set
        vav = VariantAnnotationVersion(columns_version=3)
        counted = set(vav.get_rankscore_pathogenic_prediction_funcs())
        self.assertEqual({t.rankscore_field for t in TOOLS if t.rankscore_field}, counted)

    def test_damaging_end_only_where_no_benign_cutoff(self):
        # Conservation and dbscSNV colour the damaging end - a low score isn't benign evidence
        bands = json.loads(pathogenicity_thresholds())
        for field in ("gerp_pp_rs", "phastcons_100_way_vertebrate", "dbscsnv_ada_score"):
            self.assertNotIn("benign", bands[field])

    def test_maxentscan_bands_on_magnitude(self):
        # Signed % change against the reference site - a big shift either way is worth flagging
        bands = json.loads(pathogenicity_thresholds())
        self.assertTrue(bands["maxentscan_percent_diff_ref"]["magnitude"])

    def test_cosmic_count_not_banded(self):
        # A somatic sample count isn't a pathogenicity axis (#1673)
        bands = json.loads(pathogenicity_thresholds())
        self.assertNotIn("cosmic_count", bands)
