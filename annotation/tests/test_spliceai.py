import os

from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.filters.damage_node import DamageNode
from annotation.fake_annotation import (
    create_fake_variant_annotation,
    create_fake_variants,
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import VariantAnnotation
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vep_annotation import _spliceai_label
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class TestSpliceAIMaxDsHelper(TestCase):

    def test_all_null_does_not_set_max(self):
        data = {
            "spliceai_pred_ds_ag": None,
            "spliceai_pred_ds_al": None,
            "spliceai_pred_ds_dg": None,
            "spliceai_pred_ds_dl": None,
        }
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertNotIn("spliceai_max_ds", data)

    def test_no_keys_present_does_not_set_max(self):
        data = {}
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertNotIn("spliceai_max_ds", data)

    def test_mixed_null_picks_max_of_present(self):
        data = {
            "spliceai_pred_ds_ag": 0.05,
            "spliceai_pred_ds_al": None,
            "spliceai_pred_ds_dg": 0.42,
            "spliceai_pred_ds_dl": None,
        }
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertAlmostEqual(data["spliceai_max_ds"], 0.42)

    def test_all_zero_masked_style_sets_zero(self):
        data = {
            "spliceai_pred_ds_ag": 0.0,
            "spliceai_pred_ds_al": 0.0,
            "spliceai_pred_ds_dg": 0.0,
            "spliceai_pred_ds_dl": 0.0,
        }
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertEqual(data["spliceai_max_ds"], 0.0)

    def test_raw_style_picks_highest(self):
        data = {
            "spliceai_pred_ds_ag": 0.05,
            "spliceai_pred_ds_al": 0.17,
            "spliceai_pred_ds_dg": 0.01,
            "spliceai_pred_ds_dl": 0.04,
        }
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertAlmostEqual(data["spliceai_max_ds"], 0.17)

    def test_string_values_are_coerced_to_float(self):
        data = {
            "spliceai_pred_ds_ag": "0.10",
            "spliceai_pred_ds_al": "0.30",
            "spliceai_pred_ds_dg": None,
            "spliceai_pred_ds_dl": None,
        }
        BulkVEPVCFAnnotationInserter._add_spliceai_max_ds(data)
        self.assertAlmostEqual(data["spliceai_max_ds"], 0.30)


class TestSpliceAILabel(TestCase):

    def test_raw_no_header_falls_back_to_flavour(self):
        path = "/does/not/exist/spliceai_scores.raw.snv.hg19.vcf.gz"
        self.assertEqual(_spliceai_label(path), "raw")

    def test_masked_no_header_falls_back_to_flavour(self):
        path = "/does/not/exist/spliceai_scores.masked.snv.hg38.vcf.gz"
        self.assertEqual(_spliceai_label(path), "masked")

    def test_reads_version_from_vcf_header(self):
        import tempfile

        header = (
            '##fileformat=VCFv4.1\n'
            '##INFO=<ID=SpliceAI,Number=.,Type=String,'
            'Description="SpliceAIv1.3 variant annotation. '
            'These include delta scores (DS) and delta positions (DP) for ...">\n'
            '#CHROM\tPOS\tID\tREF\tALT\n'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".raw.snv.vcf",
                                         delete=False) as fh:
            fh.write(header)
            path = fh.name
        try:
            self.assertEqual(_spliceai_label(path), "raw 1.3")
        finally:
            os.unlink(path)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestSpliceAIBackfillCommand(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.av = get_fake_annotation_version(cls.grch37)
        cls.vav = cls.av.variant_annotation_version
        create_fake_variants(cls.grch37)

    def test_backfill_populates_max_and_keeps_all_null_rows_null(self):
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:4])
        self.assertGreaterEqual(len(variants), 4)

        cases = [
            (variants[0], dict(spliceai_pred_ds_ag=0.05, spliceai_pred_ds_al=0.17,
                               spliceai_pred_ds_dg=0.01, spliceai_pred_ds_dl=0.04), 0.17),
            (variants[1], dict(spliceai_pred_ds_ag=None, spliceai_pred_ds_al=None,
                               spliceai_pred_ds_dg=0.0, spliceai_pred_ds_dl=0.0), 0.0),
            (variants[2], dict(spliceai_pred_ds_ag=0.9, spliceai_pred_ds_al=None,
                               spliceai_pred_ds_dg=None, spliceai_pred_ds_dl=None), 0.9),
            (variants[3], dict(spliceai_pred_ds_ag=None, spliceai_pred_ds_al=None,
                               spliceai_pred_ds_dg=None, spliceai_pred_ds_dl=None), None),
        ]

        vas = []
        for variant, ds_values, _ in cases:
            va = create_fake_variant_annotation(variant, self.vav)
            for k, v in ds_values.items():
                setattr(va, k, v)
            va.spliceai_max_ds = None
            va.save()
            vas.append(va)

        VariantAnnotation.backfill_spliceai_max_ds(self.vav)

        for va, (_variant, _ds, expected) in zip(vas, cases):
            va.refresh_from_db()
            if expected is None:
                self.assertIsNone(va.spliceai_max_ds)
            else:
                self.assertAlmostEqual(va.spliceai_max_ds, expected)

    def test_backfill_is_idempotent(self):
        variant = Variant.objects.filter(Variant.get_no_reference_q()).first()
        va = create_fake_variant_annotation(variant, self.vav)
        va.spliceai_pred_ds_ag = 0.42
        va.spliceai_pred_ds_al = None
        va.spliceai_pred_ds_dg = None
        va.spliceai_pred_ds_dl = None
        va.spliceai_max_ds = None
        va.save()

        VariantAnnotation.backfill_spliceai_max_ds(self.vav)
        VariantAnnotation.backfill_spliceai_max_ds(self.vav)
        va.refresh_from_db()
        self.assertAlmostEqual(va.spliceai_max_ds, 0.42)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestDamageNodeSpliceAIQ(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        from django.contrib.auth.models import User
        user = User.objects.get_or_create(username="test_TestDamageNodeSpliceAIQ")[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def test_splice_min_filters_on_spliceai_max_ds(self):
        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__spliceai_max_ds__gte", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_ag__gte", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_al__gte", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_dg__gte", q_str)
        self.assertNotIn("variantannotation__spliceai_pred_ds_dl__gte", q_str)

    def test_splice_required_with_allow_null_includes_isnull(self):
        node = DamageNode(
            analysis=self.analysis,
            splice_min=0.5,
            splice_required=True,
            splice_allow_null=True,
        )
        q_str = str(node._get_node_q())
        self.assertIn("variantannotation__spliceai_max_ds__isnull", q_str)

    def test_splice_required_without_allow_null_excludes_isnull(self):
        node = DamageNode(
            analysis=self.analysis,
            splice_min=0.5,
            splice_required=True,
            splice_allow_null=False,
        )
        q_str = str(node._get_node_q())
        self.assertNotIn("variantannotation__spliceai_max_ds__isnull", q_str)

    def test_splice_min_returns_expected_rows(self):
        from snpdb.models import Variant
        create_fake_variants(self.grch37)
        vav = self.analysis.annotation_version.variant_annotation_version

        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        self.assertGreaterEqual(len(variants), 3)

        v_high = variants[0]
        v_low = variants[1]
        v_null = variants[2]

        va_high = create_fake_variant_annotation(v_high, vav)
        va_high.spliceai_max_ds = 0.9
        va_high.save()

        va_low = create_fake_variant_annotation(v_low, vav)
        va_low.spliceai_max_ds = 0.1
        va_low.save()

        va_null = create_fake_variant_annotation(v_null, vav)
        va_null.spliceai_max_ds = None
        va_null.save()

        node = DamageNode(analysis=self.analysis, splice_min=0.5)
        q = node._get_node_q()
        matched_ids = set(
            Variant.objects.filter(q).filter(pk__in=[v_high.pk, v_low.pk, v_null.pk])
                .values_list("pk", flat=True)
        )
        self.assertIn(v_high.pk, matched_ids)
        self.assertNotIn(v_low.pk, matched_ids)
        self.assertNotIn(v_null.pk, matched_ids)
