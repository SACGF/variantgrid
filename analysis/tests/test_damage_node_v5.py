"""
DamageNode at columns_version 5 (#579): rankscores become legacy-only, and the PTC NMD filter
(ptc_nmd) joins the LoF filters. Also the gnomAD gene constraint and repeat region filters (#2600), run as SQL.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.enums import PTCNMDFilter
from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import AnnotationRun, VariantAnnotationVersion
from annotation.models.models import VariantAnnotation
from genes.models.models_gene import GeneSymbol
from genes.models.models_gnomad_gene_constraint import GnomADGeneConstraint
from library.django_utils.django_partition import temporary_db_table
from snpdb.models.models_cached_web_resource import CachedWebResource
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


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

    def test_ptc_nmd_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING)
        self.assertTrue(node.modifies_parents())
        # Categorical LoF flag, not a per-tool pathogenicity prediction
        self.assertFalse(node.has_individual_pathogenic_predictions())

    def test_ptc_nmd_emits_q(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__nmd_escape_status", q_str)
        self.assertNotIn("PREDICTED_NMD", q_str)

    def test_ptc_nmd_required_is_and_filter(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING, ptc_nmd_required=True,
                          alphamissense_score_min=0.5)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__nmd_escape_status", q_str)
        self.assertTrue(node.has_required())

    def test_ptc_nmd_escaping_or_predicted_keeps_nmd_predicted(self):
        """ #2596 - a frameshift predicted to undergo NMD is kept, not just the escaping ones """
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING_OR_PREDICTED)
        q_str = str(node._get_node_q())
        self.assertIn("PREDICTED_NMD", q_str)
        self.assertIn("ESCAPING", q_str)

    def test_tool_not_annotated_on_build_does_not_modify_parents(self):
        """ EVE is GRCh38 only, so a value left on a GRCh37 node filters nothing (#2345) """
        node = DamageNode(analysis=self.analysis, eve_score_min=0.5, eve_score_required=True)
        self.assertFalse(node.modifies_parents())
        self.assertFalse(node.has_required())

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

    def test_ptc_nmd_modifies_parents(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING)
        self.assertTrue(node.modifies_parents())


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV4IgnoresV5Fields(_ColumnsVersion5Mixin, TestCase):
    COLUMNS_VERSION = 4

    def test_v4_does_not_emit_ptc_q(self):
        node = DamageNode(analysis=self.analysis, ptc_nmd=PTCNMDFilter.ESCAPING)
        q_str = str(node._get_node_q() or "")
        self.assertNotIn("nmd_escape_status", q_str)

    @override_settings(ANNOTATION_SHOW_LEGACY_RANKSCORES=True)
    def test_v4_still_honours_legacy_rankscore_setting(self):
        node = DamageNode(analysis=self.analysis)
        self.assertTrue(node.show_rankscore_sliders())


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
@patch.object(VariantAnnotationVersion, "visible_columns", frozenset({"repeat_masker"}))
class TestDamageNodeGeneConstraintAndRepeat(_ColumnsVersion5Mixin, TestCase):
    """ #2600 - gnomAD gene constraint (PVS1 mechanism, PP2 / BP1) and PM4's non-repeat region """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        annotation_run = AnnotationRun.objects.create()
        cached_web_resource = CachedWebResource.objects.create(name="test_gnomad_gene_constraint")
        partition_table = vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        cls.variants = {}
        for i, (name, symbol, repeat_masker) in enumerate([("constrained", "CONSTRAINED1", None),
                                                          ("tolerant", "TOLERANT1", None),
                                                          ("no_constraint_row", "NOROW1", None),
                                                          ("intergenic", None, None),
                                                          ("repeat", "CONSTRAINED1", "AluY")]):
            variant = slowly_create_test_variant("1", 1000 + i * 100, "A", "T", cls.grch37)
            with temporary_db_table(VariantAnnotation, partition_table):
                VariantAnnotation.objects.create(version=vav, variant=variant, annotation_run=annotation_run,
                                                 symbol=symbol, repeat_masker=repeat_masker,
                                                 predictions_num_pathogenic=0, predictions_num_benign=0)
            cls.variants[name] = variant

        # (symbol, mane_select, LOEUF, missense Z) - TOLERANT1's non-MANE transcript is constrained, its gene isn't
        for symbol, mane_select, loeuf, mis_z in [("CONSTRAINED1", True, 0.2, 4.0),
                                                  ("TOLERANT1", True, 1.2, 0.5),
                                                  ("TOLERANT1", False, 0.1, 5.0)]:
            gene_symbol = GeneSymbol.objects.get_or_create(symbol=symbol)[0]
            GnomADGeneConstraint.objects.create(gene_symbol=gene_symbol, cached_web_resource=cached_web_resource,
                                                mane_select=mane_select, lof_oe_ci_upper=loeuf, mis_z_score=mis_z)

    def _filtered(self, node: DamageNode) -> set[str]:
        q = node._get_node_q()
        pk_names = {v.pk: name for name, v in self.variants.items()}
        qs = get_variant_queryset_for_annotation_version(self.analysis.annotation_version)
        qs = qs.filter(pk__in=pk_names).filter(q)
        return {pk_names[pk] for pk in qs.values_list("pk", flat=True)}

    def test_loeuf_uses_mane_select_row(self):
        node = DamageNode(analysis=self.analysis, gene_loeuf_max=0.6)
        self.assertEqual({"constrained", "repeat"}, self._filtered(node))

    def test_missense_z(self):
        node = DamageNode(analysis=self.analysis, gene_mis_z_min=3.09)
        self.assertEqual({"constrained", "repeat"}, self._filtered(node))

    def test_required_allow_null_keeps_genes_without_constraint(self):
        node = DamageNode(analysis=self.analysis, gene_loeuf_max=0.6, gene_loeuf_required=True,
                          gene_loeuf_allow_null=True)
        self.assertEqual({"constrained", "repeat", "no_constraint_row", "intergenic"}, self._filtered(node))
        node.gene_loeuf_allow_null = False
        self.assertEqual({"constrained", "repeat"}, self._filtered(node))

    def test_outside_repeat_region_restricts(self):
        node = DamageNode(analysis=self.analysis, gene_loeuf_max=0.6, outside_repeat_region=True)
        self.assertEqual({"constrained"}, self._filtered(node))

    def test_gene_constraint_skipped_without_import(self):
        GnomADGeneConstraint.objects.all().delete()
        node = DamageNode(analysis=self.analysis, gene_loeuf_max=0.6, gene_loeuf_required=True)
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())
        self.assertTrue([w for w in node.get_warnings() if "Gene constraint filter is NOT being applied" in w])
