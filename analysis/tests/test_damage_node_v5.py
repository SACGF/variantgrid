"""
DamageNode at columns_version 5 (#579): the version dicts gain a `5:` key with v5 as the
fallback, rankscores become legacy-only, and ptc_nmd_escaping joins the LoF filters.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild


class _ColumnsVersion5Mixin:
    """ columns_version is a derived property on AnalysisNode (reads through the VAV), so forcing
        the VAV is the only way to flip a node into its v5 code path. """
    COLUMNS_VERSION = 5

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


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV5(_ColumnsVersion5Mixin, TestCase):

    def test_raw_score_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.5)
        self.assertTrue(node.modifies_parents())
        self.assertTrue(node.has_individual_pathogenic_predictions())

    def test_pred_field_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, alphamissense_pred="P")
        self.assertTrue(node.modifies_parents())

    def test_ptc_nmd_escaping_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd_escaping=True)
        self.assertTrue(node.modifies_parents())
        # Categorical LoF flag, not a per-tool pathogenicity prediction
        self.assertFalse(node.has_individual_pathogenic_predictions())

    def test_ptc_nmd_escaping_emits_q(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd_escaping=True)
        self.assertIn("variantannotation__nmd_escape_status", str(node._get_node_q()))

    def test_ptc_nmd_escaping_required_is_and_filter(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd_escaping=True, ptc_nmd_escaping_required=True,
                          alphamissense_score_min=0.5)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__nmd_escape_status", q_str)
        self.assertTrue(node.has_required())

    def test_empty_node_does_not_modify_parents(self):
        node = DamageNode(analysis=self.analysis)
        self.assertFalse(node.modifies_parents())
        self.assertFalse(node.has_required())

    def test_rankscore_set_still_filters_and_shows_sliders(self):
        """ Legacy rankscore filters keep applying at v5, so the node keeps reporting them. """
        node = DamageNode(analysis=self.analysis, revel_rankscore_min=0.5, revel_rankscore_required=True)
        self.assertTrue(node.show_rankscore_sliders())
        self.assertTrue(node.modifies_parents())
        self.assertTrue(node.has_required())
        self.assertIn("variantannotation__revel_rankscore__gte", str(node._get_node_q()))

    def test_no_rankscore_hides_sliders(self):
        node = DamageNode(analysis=self.analysis)
        self.assertFalse(node.show_rankscore_sliders())

    @override_settings(ANNOTATION_SHOW_LEGACY_RANKSCORES=True)
    def test_legacy_rankscore_setting_does_not_reveal_sliders_at_v5(self):
        """ ANNOTATION_SHOW_LEGACY_RANKSCORES keeps its meaning at v4 and stops offering
            rankscores once a deployment reaches v5. """
        node = DamageNode(analysis=self.analysis)
        self.assertFalse(node.show_rankscore_sliders())


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV6InheritsV5(_ColumnsVersion5Mixin, TestCase):
    """ The version dicts fall back to the v5 set, so a future columns_version isn't inert. """
    COLUMNS_VERSION = 6

    def test_raw_score_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.5)
        self.assertTrue(node.modifies_parents())

    def test_ptc_nmd_escaping_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd_escaping=True)
        self.assertTrue(node.modifies_parents())


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV4IgnoresV5Fields(_ColumnsVersion5Mixin, TestCase):
    COLUMNS_VERSION = 4

    def test_v4_does_not_emit_ptc_q(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd_escaping=True)
        q_str = str(node._get_node_q() or "")
        self.assertNotIn("nmd_escape_status", q_str)

    @override_settings(ANNOTATION_SHOW_LEGACY_RANKSCORES=True)
    def test_v4_still_honours_legacy_rankscore_setting(self):
        node = DamageNode(analysis=self.analysis)
        self.assertTrue(node.show_rankscore_sliders())
