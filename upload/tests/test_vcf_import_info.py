from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from library.utils import sha256sum_str
from snpdb.models import Sequence, GenomeBuild
from snpdb.models.models_enums import ImportSource
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from upload.models import (
    UploadedFile, UploadPipeline, UploadStep, UploadedFileTypes,
    SimpleVCFImportInfo, ModifiedImportedVariants, ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
)


def _make_upload_step(user, step_name="test step"):
    uploaded_file = UploadedFile.objects.create(
        user=user, name="test.vcf", path="/tmp/test.vcf",
        file_type=UploadedFileTypes.VCF, import_source=ImportSource.COMMAND_LINE,
    )
    pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
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
