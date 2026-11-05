"""
DamageNode filters over the plugins added to annotation versions in use (#2561): AVI, AlphaGenome splicing,
SpliceVault and UTRAnnotator - skipped until the version is backfilled, or where it withholds the column (#2571),
AlphaGenome splicing joining the splicing OR.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.enums import SpliceVaultOutOfFrame, UTRHighImpact
from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import create_annotated_variants, node_keeps
from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild

PLUGIN_COLUMNS = frozenset({"avi_phred", "alphagenome_splicing", "splicevault_out_of_frame_events",
                            "utr_annotator_consequence", "utr_annotator_annotation"})
PLUGIN_FILTERS = {"avi_phred_min": 20, "alphagenome_splicing_min": 0.5,
                  "splicevault_out_of_frame": SpliceVaultOutOfFrame.ALL, "utr_high_impact": UTRHighImpact.ANY_FRAME}

UAUG_GAIN = "5_prime_UTR_premature_start_codon_gain_variant"
_UAUG_GAIN_RECORD = "CapDistanceToStart=111:DistanceToCDS=4:DistanceToStop=126:KozakContext=GTCATGG:KozakStrength={}"
UTR_HIGH_IMPACT = {"utr_annotator_consequence": UAUG_GAIN,
                   "utr_annotator_annotation": _UAUG_GAIN_RECORD.format("Strong:type=OutOfFrame_oORF")}
UTR_WEAK_KOZAK = {"utr_annotator_consequence": UAUG_GAIN,
                  "utr_annotator_annotation": _UAUG_GAIN_RECORD.format("Weak:type=OutOfFrame_oORF")}
SPLICEVAULT_TOP_OUT_OF_FRAME = {"splicevault_top_events": "Top1:ES:2:45%:Frameshift&Top2:CD:-14:10%:Inframe",
                                "splicevault_out_of_frame_events": "Frameshift:1/2"}
SPLICEVAULT_TOP_IN_FRAME = {"splicevault_top_events": "Top1:ES:2:45%:Inframe&Top2:CD:-14:10%:Frameshift",
                            "splicevault_out_of_frame_events": "Frameshift:1/2"}
# dbscSNV is in the splicing OR too, so it's set (and failing) on every splicing row
SPLICE_SCORES_FAIL = {"spliceai_max_ds": 0.1, "dbscsnv_ada_score": 0.1, "dbscsnv_rf_score": 0.1}

ROWS = {
    # Each plugin's PLUGIN_FILTERS threshold - one passing row apiece, and one failing them all
    "avi": {"avi_phred": 25},
    "splicevault_all": {"splicevault_out_of_frame_events": "Frameshift:2/2"},
    "utr": UTR_HIGH_IMPACT,
    "fails_every_plugin": {"avi_phred": 10, "alphagenome_splicing": 0.2,
                           "splicevault_out_of_frame_events": "Frameshift:1/4", **UTR_WEAK_KOZAK},
    # Splicing: SpliceAI 0.5 OR AlphaGenome 1.0
    "spliceai_only": {**SPLICE_SCORES_FAIL, "spliceai_max_ds": 0.8, "alphagenome_splicing": 0.2},
    "alphagenome_only": {**SPLICE_SCORES_FAIL, "alphagenome_splicing": 1.2},
    "null_alphagenome": SPLICE_SCORES_FAIL,
    "fails_both_splicing": {**SPLICE_SCORES_FAIL, "alphagenome_splicing": 0.5},
    # SpliceVault top event out of frame AND out-of-frame high-impact UTR
    "splicevault_and_utr": {**SPLICEVAULT_TOP_OUT_OF_FRAME, **UTR_HIGH_IMPACT},
    "null_splicevault": UTR_HIGH_IMPACT,
    "null_utr": SPLICEVAULT_TOP_OUT_OF_FRAME,
    "splicevault_in_frame": {**SPLICEVAULT_TOP_IN_FRAME, **UTR_HIGH_IMPACT},
}


# The fake version has no AlphaGenome pin, so it would withhold that column: the tests say which columns are offered
@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
@patch.object(VariantAnnotationVersion, "visible_columns", PLUGIN_COLUMNS)
@patch.object(VariantAnnotationVersion, "withheld_licence_restricted_columns", frozenset())
class DamageNodePluginFiltersTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        vav = VariantAnnotationVersion.objects.filter(genome_build=grch37).first()
        vav.columns_version = 5
        vav.save()
        user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.analysis = Analysis(genome_build=grch37)
        cls.analysis.set_defaults_and_save(user)

        cls.variants = create_annotated_variants(vav, ROWS)

    def _filtered(self, node: DamageNode, among: set[str]) -> set[str]:
        return node_keeps(node, self.variants, among)

    def _set_backfilled(self, backfilled: bool):
        vav = self.analysis.annotation_version.variant_annotation_version
        vav.backfilled_utr_avi_splicevault = backfilled
        vav.save()

    def test_skipped_with_warning_until_backfilled(self):
        """ Null on not-yet-backfilled variants means "not imported", so filtering would drop them all """
        self._set_backfilled(False)
        node = DamageNode(analysis=self.analysis, **PLUGIN_FILTERS)
        self.assertTrue(node.modifies_parents())
        self.assertIsNone(node._get_node_q())
        not_applied = [w for w in node.get_warnings() if "NOT being applied" in w]
        self.assertEqual(4, len(not_applied))

    def test_applied_once_backfilled(self):
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, **PLUGIN_FILTERS)
        among = {"avi", "alphagenome_only", "splicevault_all", "utr", "fails_every_plugin"}
        self.assertEqual({"avi", "alphagenome_only", "splicevault_all", "utr"}, self._filtered(node, among))
        self.assertFalse([w for w in node.get_warnings() if "NOT being applied" in w])

    def test_withheld_non_commercial_filter_not_applied(self):
        """ #2571: a saved AlphaGenome threshold on a version that withholds it (licence rules out the deployment's use) mustn't still filter """
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, alphagenome_splicing_min=0.5, avi_phred_min=20)
        withheld = frozenset({"alphagenome_splicing"})
        with (patch.object(VariantAnnotationVersion, "visible_columns", PLUGIN_COLUMNS - withheld),
              patch.object(VariantAnnotationVersion, "withheld_licence_restricted_columns", withheld)):
            q_str = str(node._get_node_q())
            warnings = node.get_warnings()
        self.assertNotIn("alphagenome_splicing", q_str)
        self.assertIn("avi_phred__gte", q_str)
        self.assertEqual(1, len([w for w in warnings if "AlphaGenome splicing filter is NOT being applied" in w]))

    def test_alphagenome_splicing_is_a_splicing_or_term(self):
        """ Required splicing keeps a variant passing SpliceAI or AlphaGenome, or with no AlphaGenome score """
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, splice_min=0.5, alphagenome_splicing_min=1.0,
                          splice_required=True, splice_allow_null=True)
        among = {"spliceai_only", "alphagenome_only", "null_alphagenome", "fails_both_splicing"}
        self.assertEqual({"spliceai_only", "alphagenome_only", "null_alphagenome"}, self._filtered(node, among))

    def test_alphagenome_splicing_without_splice_min(self):
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, alphagenome_splicing_min=1.0)
        among = {"spliceai_only", "alphagenome_only", "null_alphagenome", "fails_both_splicing"}
        self.assertEqual({"alphagenome_only"}, self._filtered(node, among))

    def test_required_allow_null(self):
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis,
                          splicevault_out_of_frame=SpliceVaultOutOfFrame.TOP, splicevault_out_of_frame_required=True,
                          splicevault_out_of_frame_allow_null=False,
                          utr_high_impact=UTRHighImpact.OUT_OF_FRAME, utr_high_impact_required=True,
                          utr_high_impact_allow_null=False)
        among = {"splicevault_and_utr", "null_splicevault", "null_utr", "splicevault_in_frame"}
        self.assertEqual({"splicevault_and_utr"}, self._filtered(node, among))

        node.splicevault_out_of_frame_allow_null = True
        node.utr_high_impact_allow_null = True
        self.assertEqual({"splicevault_and_utr", "null_splicevault", "null_utr"}, self._filtered(node, among))
