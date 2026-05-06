from django.test import TestCase, override_settings

from annotation.fake_annotation import (
    create_fake_variant_annotation,
    create_fake_variants,
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import VariantAnnotation
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


V4_FIELDS = ["af_1kg", "af_uk10k", "gnomad_af", "gnomad_popmax_af", "gnomad_fafmax_faf95_max"]
V3_FIELDS = ["af_1kg", "af_uk10k", "gnomad_af", "gnomad_popmax_af"]


class _StubInserter:
    def __init__(self, fields):
        self._max_af_field_set = fields

    _add_max_af = BulkVEPVCFAnnotationInserter._add_max_af


class TestMaxAfHelper(TestCase):

    def test_v4_all_populated_picks_max(self):
        data = {
            "af_1kg": 0.01,
            "af_uk10k": 0.02,
            "gnomad_af": 0.05,
            "gnomad_popmax_af": 0.07,
            "gnomad_fafmax_faf95_max": 0.06,
        }
        _StubInserter(V4_FIELDS)._add_max_af(data)
        self.assertAlmostEqual(data["max_af"], 0.07)

    def test_v4_one_field_null_leaves_max_unset(self):
        for missing in V4_FIELDS:
            data = {f: 0.01 for f in V4_FIELDS}
            data[missing] = None
            _StubInserter(V4_FIELDS)._add_max_af(data)
            self.assertNotIn("max_af", data, f"max_af should be unset when {missing} is NULL")

    def test_v4_all_null_leaves_max_unset(self):
        data = {f: None for f in V4_FIELDS}
        _StubInserter(V4_FIELDS)._add_max_af(data)
        self.assertNotIn("max_af", data)

    def test_v4_one_field_missing_key_leaves_max_unset(self):
        data = {f: 0.01 for f in V4_FIELDS}
        del data["gnomad_fafmax_faf95_max"]
        _StubInserter(V4_FIELDS)._add_max_af(data)
        self.assertNotIn("max_af", data)

    def test_v3_all_populated_picks_max(self):
        data = {
            "af_1kg": 0.01,
            "af_uk10k": 0.02,
            "gnomad_af": 0.05,
            "gnomad_popmax_af": 0.07,
        }
        _StubInserter(V3_FIELDS)._add_max_af(data)
        self.assertAlmostEqual(data["max_af"], 0.07)

    def test_v3_popmax_null_leaves_max_unset(self):
        data = {
            "af_1kg": 0.01,
            "af_uk10k": 0.02,
            "gnomad_af": 0.05,
            "gnomad_popmax_af": None,
        }
        _StubInserter(V3_FIELDS)._add_max_af(data)
        self.assertNotIn("max_af", data)

    def test_string_values_are_coerced_to_float(self):
        data = {
            "af_1kg": "0.01",
            "af_uk10k": "0.02",
            "gnomad_af": "0.05",
            "gnomad_popmax_af": "0.07",
        }
        _StubInserter(V3_FIELDS)._add_max_af(data)
        self.assertAlmostEqual(data["max_af"], 0.07)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class TestMaxAfBackfillCommand(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.av = get_fake_annotation_version(cls.grch37)
        cls.vav = cls.av.variant_annotation_version
        create_fake_variants(cls.grch37)

    def test_backfill_populates_full_info_skips_partial(self):
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        self.assertGreaterEqual(len(variants), 3)

        full_values = dict(af_1kg=0.01, af_uk10k=0.02, gnomad_af=0.05, gnomad_popmax_af=0.07)
        partial_values = dict(af_1kg=0.01, af_uk10k=0.02, gnomad_af=0.05, gnomad_popmax_af=None)
        all_null_values = dict(af_1kg=None, af_uk10k=None, gnomad_af=None, gnomad_popmax_af=None)

        cases = [
            (variants[0], full_values, 0.07),
            (variants[1], partial_values, None),
            (variants[2], all_null_values, None),
        ]

        vas = []
        for variant, af_values, _ in cases:
            va = create_fake_variant_annotation(variant, self.vav)
            for k, v in af_values.items():
                setattr(va, k, v)
            va.max_af = None
            va.save()
            vas.append(va)

        VariantAnnotation.backfill_max_af(self.vav)

        for va, (_variant, _af, expected) in zip(vas, cases):
            va.refresh_from_db()
            if expected is None:
                self.assertIsNone(va.max_af)
            else:
                self.assertAlmostEqual(va.max_af, expected)

    def test_backfill_is_idempotent(self):
        variant = Variant.objects.filter(Variant.get_no_reference_q()).first()
        va = create_fake_variant_annotation(variant, self.vav)
        va.af_1kg = 0.01
        va.af_uk10k = 0.02
        va.gnomad_af = 0.05
        va.gnomad_popmax_af = 0.07
        va.max_af = None
        va.save()

        first = VariantAnnotation.backfill_max_af(self.vav)
        second = VariantAnnotation.backfill_max_af(self.vav)
        self.assertGreaterEqual(first, 1)
        self.assertEqual(second, 0)
        va.refresh_from_db()
        self.assertAlmostEqual(va.max_af, 0.07)
