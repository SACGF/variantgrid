import random

from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from library.utils import sha256sum_str
from snpdb.models import GenomeBuild, Sequence
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_variant import Variant, VariantCoordinate
from snpdb.signals.variant_search import get_results_from_variant_coordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    ModifiedImportedVariants,
    SimpleVCFImportInfo,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
)


def _make_upload_step(user, step_name="test step"):
    file_upload = FileUpload.objects.create(
        user=user, name="test.vcf", path="/tmp/test.vcf",
        file_type=UploadedFileTypes.VCF, import_source=ImportSource.COMMAND_LINE,
    )
    pipeline = UploadPipeline.objects.create(file_upload=file_upload)
    return UploadStep.objects.create(upload_pipeline=pipeline, name=step_name, sort_order=0)


class TestSimpleVCFImportInfoAddMessageCount(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.create_user(username="info_test_user", password="x")
        cls.upload_step = _make_upload_step(user)

    def test_accumulates_count_on_second_call(self):
        """Two serial calls for the same message must merge into one row, not create two."""
        msg = "accumulate_test_unique_xyz"
        SimpleVCFImportInfo.add_message_count(5, msg, self.upload_step,
                                              type=SimpleVCFImportInfo.SVLEN_MODIFIED, has_more_details=True)
        SimpleVCFImportInfo.add_message_count(3, msg, self.upload_step,
                                              type=SimpleVCFImportInfo.SVLEN_MODIFIED, has_more_details=True)
        qs = SimpleVCFImportInfo.objects.filter(message_string=msg)
        self.assertEqual(qs.count(), 1, "Two calls for same message should produce exactly 1 row")
        self.assertEqual(qs.first().count, 8, "Count should be 5 + 3 = 8")


class TestModifiedImportedVariantsMessage(TestCase):
    """Tests for ModifiedImportedVariants.message — focuses on the Postgres regexp_replace
    deduplication logic, which is the only non-trivial part of the property."""

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_sha256_hash=sha256sum_str(base))
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        user = User.objects.create_user(username="mivs_test_user", password="x")
        upload_step = _make_upload_step(user, step_name="Normalise variants")
        cls.mivs = ModifiedImportedVariants.objects.create(upload_step=upload_step)
        cls.variant1 = slowly_create_test_variant("1", 100, "A", "T", cls.grch37)
        cls.variant2 = slowly_create_test_variant("1", 200, "C", "G", cls.grch37)

    def test_message_multiallelic_deduplication(self):
        """Two MIVs from the same multi-allelic row (different alt indices) count as 1, not 2.

        The message property strips the trailing |N index via Postgres regexp_replace and then
        runs distinct(). Records "1|100|A|C,T|1" and "1|100|A|C,T|2" share the stripped form
        "1|100|A|C,T" and should appear as one multi-allelic event in the report.
        """
        ModifiedImportedVariant.objects.create(
            import_info=self.mivs, variant=self.variant1,
            operation=ModifiedImportedVariantOperation.NORMALIZATION,
            old_multiallelic="1|100|A|C,T|1", old_variant=None, old_variant_formatted="1:100:A/C",
        )
        ModifiedImportedVariant.objects.create(
            import_info=self.mivs, variant=self.variant2,
            operation=ModifiedImportedVariantOperation.NORMALIZATION,
            old_multiallelic="1|100|A|C,T|2", old_variant=None, old_variant_formatted="1:100:A/T",
        )
        msg = self.mivs.message
        self.assertIn("1 multi-allelic split", msg, "Two alts from same row should count as 1 event")
        self.assertNotIn("2 multi-allelic split", msg)


class TestModifiedImportedVariantLongValueLookup(TestCase):
    """ Long indels make old values bigger than a btree entry can hold - lookups go via an indexed prefix,
        so values sharing that prefix must still be told apart by the full value """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_sha256_hash=sha256sum_str(base))
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        user = User.objects.create_user(username="miv_long_test_user", password="x")
        mivs = ModifiedImportedVariants.objects.create(upload_step=_make_upload_step(user))
        cls.variant1 = slowly_create_test_variant("1", 100, "A", "T", cls.grch37)
        cls.variant2 = slowly_create_test_variant("1", 200, "C", "G", cls.grch37)

        # Random so it doesn't compress under the index entry size limit
        cls.long_ref = "".join(random.Random(2092).choices("GATC", k=40000))
        old_multiallelic = f"1|100|{cls.long_ref}|C,T|"
        for i, (variant, alt) in enumerate([(cls.variant1, "C"), (cls.variant2, "T")], start=1):
            ModifiedImportedVariant.objects.create(
                import_info=mivs, variant=variant, operation=ModifiedImportedVariantOperation.NORMALIZATION,
                old_multiallelic=f"{old_multiallelic}{i}",
                old_variant_formatted=f"1:100:{cls.long_ref}/{alt}",
            )

    def test_exact_lookup_uses_full_value(self):
        vc = VariantCoordinate(chrom="1", position=100, ref=self.long_ref, alt="T")
        variants = ModifiedImportedVariant.get_variants_for_unnormalized_variant(Variant.objects.all(), vc)
        self.assertEqual(list(variants), [self.variant2])

    def test_any_alt_lookup(self):
        vc = VariantCoordinate(chrom="1", position=100, ref=self.long_ref, alt="")
        variants = ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(Variant.objects.all(), vc)
        self.assertEqual(set(variants), {self.variant1, self.variant2})

    def test_search_fallback_stays_inside_searched_variants(self):
        hidden_variant2_qs = Variant.objects.exclude(pk=self.variant2.pk)
        vc = VariantCoordinate(chrom="1", position=100, ref=self.long_ref, alt="T")
        self.assertFalse(get_results_from_variant_coordinate(self.grch37, hidden_variant2_qs, vc).exists())
        vc_any_alt = VariantCoordinate(chrom="1", position=100, ref=self.long_ref, alt="")
        any_alt = get_results_from_variant_coordinate(self.grch37, hidden_variant2_qs, vc_any_alt, any_alt=True)
        self.assertEqual(list(any_alt), [self.variant1])


class TestModifiedImportedVariantNormalised(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_sha256_hash=sha256sum_str(base))
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)

        user = User.objects.create_user(username="miv_normalised_test_user", password="x")
        mivs = ModifiedImportedVariants.objects.create(upload_step=_make_upload_step(user))
        variant = slowly_create_test_variant("1", 100, "A", "T", grch37)
        # Imports before #2126 left old_variant null when bcftools trimmed an indel in place
        cls.trimmed = ModifiedImportedVariant.objects.create(
            import_info=mivs, variant=variant, operation=ModifiedImportedVariantOperation.NORMALIZATION,
            old_variant_formatted="1:100:AC/TC")
        cls.split_only = ModifiedImportedVariant.objects.create(
            import_info=mivs, variant=variant, operation=ModifiedImportedVariantOperation.NORMALIZATION,
            old_multiallelic="1|100|A|G,T|2", old_variant_formatted="1:100:A/T")

    def test_q_normalised_compares_with_variant(self):
        qs = ModifiedImportedVariant.objects.filter(ModifiedImportedVariant.q_normalised())
        self.assertEqual(list(qs), [self.trimmed])

    def test_variant_coordinate_from_old_variant_alt_contig(self):
        vc = ModifiedImportedVariant.get_variant_coordinate_from_old_variant("HLA-A*01:01:01:01:100:AC/A")
        self.assertEqual(vc.as_tuple, ("HLA-A*01:01:01:01", 100, "AC", "A", None))
