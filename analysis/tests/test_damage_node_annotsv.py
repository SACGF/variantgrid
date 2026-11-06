"""
Effect node AnnotSV filters (#2701): SVs and small variants go through the one node, so every AnnotSV filter has to
leave small variants (no AnnotSV annotation) to the node's other filters, and do nothing on a build without AnnotSV.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.node_presets import damage_node_presets
from analysis.tests.utils import create_annotated_variants, node_keeps
from annotation.enums import Pathogenicity, PathogenicityImpact, VariantAnnotationPipelineType
from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.models import VariantAnnotationVersion
from annotation.models.models_annotation_run import AnnotationPipelineVersion
from snpdb.models.models_genome import GenomeBuild

HIGH = PathogenicityImpact.HIGH
MODIFIER = PathogenicityImpact.MODIFIER


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeAnnotSV(TestCase):
    ROWS = {
        "small_high": {"impact": HIGH},
        "small_modifier": {"impact": MODIFIER},
        # A partial exonic deletion VEP calls MODIFIER that AnnotSV ranks pathogenic
        "sv_pathogenic": {"impact": MODIFIER, "annotsv_acmg_class": Pathogenicity.PATHOGENIC,
                          "annotsv_acmg_score": 1.2},
        "sv_uncertain": {"impact": MODIFIER, "annotsv_acmg_class": Pathogenicity.UNCERTAIN,
                         "annotsv_acmg_score": 0.3},
        "sv_likely_pathogenic_segdup": {"impact": MODIFIER, "annotsv_acmg_class": Pathogenicity.LIKELY_PATHOGENIC,
                                        "annotsv_acmg_score": 0.95, "annotsv_segdup_right": "chr1:1-1000"},
    }

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        vav = VariantAnnotationVersion.objects.filter(genome_build=cls.grch37).first()
        user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)
        cls.variants = create_annotated_variants(vav, cls.ROWS)
        cls.annotsv = AnnotationPipelineVersion.objects.create(
            pipeline_type=VariantAnnotationPipelineType.ANNOTSV, genome_build=cls.grch37, code_version="3.5.8",
            status=AnnotationPipelineVersion.Status.ACTIVE)

    def _filtered(self, node: DamageNode) -> set[str]:
        return node_keeps(node, self.variants)

    def test_required_allow_null_filters_svs_only(self):
        node = DamageNode(analysis=self.analysis, annotsv_acmg_class_min=Pathogenicity.LIKELY_PATHOGENIC,
                          annotsv_acmg_class_required=True, annotsv_acmg_class_allow_null=True)
        self.assertEqual({"small_high", "small_modifier", "sv_pathogenic", "sv_likely_pathogenic_segdup"},
                         self._filtered(node))

    def test_or_pool_with_impact(self):
        node = DamageNode(analysis=self.analysis, impact_min=HIGH, annotsv_acmg_score_min=0.9)
        self.assertEqual({"small_high", "sv_pathogenic", "sv_likely_pathogenic_segdup"}, self._filtered(node))

    def test_exclude_sv_artefact_regions(self):
        node = DamageNode(analysis=self.analysis, exclude_sv_artefact_regions=True)
        self.assertEqual(set(self.ROWS) - {"sv_likely_pathogenic_segdup"}, self._filtered(node))

    def test_pvs1_preset_keeps_annotsv_pathogenic(self):
        pvs1 = next(p for p in damage_node_presets(DamageNode(analysis=self.analysis)) if p.key == "pvs1")
        node = DamageNode(analysis=self.analysis, **pvs1.values)
        self.assertEqual({"small_high", "sv_pathogenic", "sv_likely_pathogenic_segdup"}, self._filtered(node))

    def test_skipped_without_active_pipeline(self):
        self.annotsv.status = AnnotationPipelineVersion.Status.HISTORICAL
        self.annotsv.save()
        node = DamageNode(analysis=self.analysis, annotsv_acmg_class_min=Pathogenicity.LIKELY_PATHOGENIC,
                          annotsv_acmg_class_required=True, exclude_sv_artefact_regions=True)
        self.assertFalse(node.modifies_parents())
        self.assertEqual([], node.get_extra_grid_columns())
        self.assertTrue([w for w in node.get_warnings() if "AnnotSV filters are NOT being applied" in w])
        pvs1 = next(p for p in damage_node_presets(DamageNode(analysis=self.analysis)) if p.key == "pvs1")
        self.assertTrue(pvs1.values["impact_required"])

    def test_no_pipeline_lookup_without_annotsv_settings(self):
        node = DamageNode(analysis=self.analysis, impact_min=HIGH)
        with patch.object(AnnotationPipelineVersion, "get_active") as get_active:
            node.modifies_parents()
            node.has_required()
            node._get_node_q()
            node.get_warnings()
        get_active.assert_not_called()
