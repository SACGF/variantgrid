"""
DamageNode filters over the plugins added to annotation versions in use (#2561): AVI, AlphaGenome splicing,
SpliceVault and UTRAnnotator - skipped until the version is backfilled, AlphaGenome splicing joining the splicing OR.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.db.models import Q
from django.test import TestCase, override_settings

from analysis.enums import SpliceVaultOutOfFrame, UTRHighImpact
from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.effect_filters import splicevault_out_of_frame_q, utr_high_impact_q
from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild

PLUGIN_COLUMNS = frozenset({"avi_phred", "alphagenome_splicing", "splicevault_out_of_frame_events",
                            "utr_annotator_consequence", "utr_annotator_annotation"})
PLUGIN_FILTERS = {"avi_phred_min": 20, "alphagenome_splicing_min": 0.5,
                  "splicevault_out_of_frame": SpliceVaultOutOfFrame.ALL, "utr_high_impact": UTRHighImpact.ANY_FRAME}


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
@patch.object(VariantAnnotationVersion, "visible_columns", PLUGIN_COLUMNS)
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
        q_str = str(node._get_node_q())
        for lookup in ["avi_phred__gte", "alphagenome_splicing__gte", "splicevault_out_of_frame_events__regex",
                       "utr_annotator_annotation__iregex"]:
            self.assertIn(lookup, q_str)
        self.assertFalse([w for w in node.get_warnings() if "NOT being applied" in w])

    def test_alphagenome_splicing_is_a_splicing_or_term(self):
        """ Required splicing keeps a variant passing SpliceAI or AlphaGenome, or with no AlphaGenome score """
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, splice_min=0.5, alphagenome_splicing_min=1.0,
                          splice_required=True, splice_allow_null=True)
        q = node._get_node_q()
        self.assertEqual(Q.OR, q.connector)
        self.assertIn(("variantannotation__spliceai_max_ds__gte", 0.5), q.children)
        self.assertIn(("variantannotation__alphagenome_splicing__gte", 1.0), q.children)
        self.assertIn(("variantannotation__alphagenome_splicing__isnull", True), q.children)

    def test_alphagenome_splicing_without_splice_min(self):
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis, alphagenome_splicing_min=1.0)
        self.assertEqual(Q(variantannotation__alphagenome_splicing__gte=1.0), node._get_node_q())

    def test_required_allow_null(self):
        self._set_backfilled(True)
        node = DamageNode(analysis=self.analysis,
                          splicevault_out_of_frame=SpliceVaultOutOfFrame.TOP, splicevault_out_of_frame_required=True,
                          splicevault_out_of_frame_allow_null=False,
                          utr_high_impact=UTRHighImpact.OUT_OF_FRAME, utr_high_impact_required=True, utr_high_impact_allow_null=False)
        q_no_nulls = node._get_node_q()
        self.assertEqual(splicevault_out_of_frame_q(top_event_only=True) & utr_high_impact_q(out_of_frame_only=True), q_no_nulls)

        node.splicevault_out_of_frame_allow_null = True
        node.utr_high_impact_allow_null = True
        q = node._get_node_q()
        expected = ((splicevault_out_of_frame_q(top_event_only=True)
                     | Q(variantannotation__splicevault_out_of_frame_events__isnull=True))
                    & (utr_high_impact_q(out_of_frame_only=True) | Q(variantannotation__utr_annotator_consequence__isnull=True)))
        self.assertEqual(expected, q)
