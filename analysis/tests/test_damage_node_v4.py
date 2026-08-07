from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.fake_annotation import (
    create_fake_variant_annotation,
    create_fake_variants,
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import VariantAnnotationVersion
from annotation.pathogenicity_predictions import TOOLS
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV4Q(TestCase):
    """ Verifies DamageNode raw-score and pred-field filtering at columns_version 4. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        # Force VAV onto v4 so DamageNode v4 branches activate. columns_version is a
        # derived property on AnalysisNode (reads through the VAV), so this is the
        # only way to flip the node into its v4 code path.
        vav = VariantAnnotationVersion.objects.filter(genome_build=cls.grch37).first()
        vav.columns_version = 4
        vav.save()
        user = User.objects.get_or_create(username="test_TestDamageNodeV4Q")[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def test_raw_score_min_emits_per_tool_q(self):
        # Pre-existing falsy-zero idiom in DamageNode means raw_min=0 is dropped;
        # pick a small positive value (or the BayesDel midpoint) per tool.
        kwargs = {}
        expected_fields = []
        for tool in TOOLS:
            field = f"{tool.raw_field}_min" if tool.raw_field else None
            if not field or not hasattr(DamageNode, field):
                continue
            kwargs[field] = max(tool.raw_min + tool.raw_step, 0.01)
            expected_fields.append(tool.raw_field)
        node = DamageNode(analysis=self.analysis, **kwargs)
        q_str = str(node._get_node_q())
        for raw_field in expected_fields:
            self.assertIn(f"variantannotation__{raw_field}__gte", q_str)

    def test_pred_field_emits_per_tool_q(self):
        node = DamageNode(
            analysis=self.analysis,
            alphamissense_pred="P",
            clinpred_pred="D",
            metarnn_pred="D",
            primateai_pred="D",
        )
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__alphamissense_pred", q_str)
        self.assertIn("variantannotation__clinpred_pred", q_str)
        self.assertIn("variantannotation__metarnn_pred", q_str)
        self.assertIn("variantannotation__primateai_pred", q_str)

    def test_raw_score_required_with_allow_null_includes_isnull(self):
        node = DamageNode(
            analysis=self.analysis,
            revel_rankscore_min=0.5,  # ensure required logic path executed
            alphamissense_score_min=0.5,
            alphamissense_score_required=True,
            alphamissense_score_allow_null=True,
        )
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__alphamissense_score__isnull", q_str)

    def test_pred_filter_returns_expected_rows(self):
        create_fake_variants(self.grch37)
        vav = self.analysis.annotation_version.variant_annotation_version
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        self.assertGreaterEqual(len(variants), 3)

        v_path, v_benign, v_null = variants

        va_path = create_fake_variant_annotation(v_path, vav)
        va_path.alphamissense_pred = "P"
        va_path.save()
        va_benign = create_fake_variant_annotation(v_benign, vav)
        va_benign.alphamissense_pred = "B"
        va_benign.save()
        create_fake_variant_annotation(v_null, vav)  # alphamissense_pred stays null

        node = DamageNode(analysis=self.analysis, alphamissense_pred="P")
        q = node._get_node_q()
        matched = set(Variant.objects.filter(q)
                                     .filter(pk__in=[v.pk for v in variants])
                                     .values_list("pk", flat=True))
        self.assertEqual(matched, {v_path.pk})

    def test_raw_score_filter_returns_expected_rows(self):
        create_fake_variants(self.grch37)
        vav = self.analysis.annotation_version.variant_annotation_version
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        self.assertGreaterEqual(len(variants), 3)

        v_path, v_benign, v_null = variants

        va_path = create_fake_variant_annotation(v_path, vav)
        va_path.alphamissense_score = 0.9  # above 0.170 PP3-supporting
        va_path.save()
        va_benign = create_fake_variant_annotation(v_benign, vav)
        va_benign.alphamissense_score = 0.05  # below threshold
        va_benign.save()
        create_fake_variant_annotation(v_null, vav)  # alphamissense_score stays null

        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.170)
        q = node._get_node_q()
        matched = set(Variant.objects.filter(q)
                                     .filter(pk__in=[v.pk for v in variants])
                                     .values_list("pk", flat=True))
        self.assertEqual(matched, {v_path.pk})


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeV3IgnoresV4Fields(TestCase):
    """ A v3 (or earlier) DamageNode must ignore v4-only model fields, even if set. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        vav = VariantAnnotationVersion.objects.filter(genome_build=cls.grch37).first()
        vav.columns_version = 3
        vav.save()
        user = User.objects.get_or_create(username="test_TestDamageNodeV3IgnoresV4Fields")[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def test_v3_does_not_emit_v4_q(self):
        node = DamageNode(analysis=self.analysis, alphamissense_score_min=0.5)
        q_str = str(node._get_node_q() or "")
        self.assertNotIn("variantannotation__alphamissense_score__gte", q_str)
