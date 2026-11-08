""" Manual Variant Entry "create as sample" (#2719): the entries become a VCF with one sample, linked to
    the collection and kept out of the global zygosity counts. The default path stays variants only """
import os
import tempfile
from unittest.mock import patch

import cyvcf2
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.fake_data import get_fake_annotation_version
from annotation.manual_variant_entry import create_manual_variants
from annotation.models.models_manual_variant import (
    CreatedManualVariant,
    ManualVariantEntry,
    ManualVariantEntryCollection,
)
from snpdb.enums import ImportSource
from snpdb.forms import ManualVariantEntryForm
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.vcf_utils import write_vcf_from_variant_coordinates
from upload.enums import UploadedFileTypes
from upload.models import (
    FileUpload,
    UploadedManualVariantEntryCollection,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.vcf.genotype_vcf_tasks import ImportCreateVCFModelForManualVariantEntryTask
from upload.vcf.bulk_manual_variant_entry_linking_vcf_processor import (
    BulkManualVariantEntryLinkingVCFProcessor,
)

VARIANTS_TEXT = "1:169519049 T>C\n13:32914438 T>TA"


@override_settings(UPLOAD_ENABLED=True, VARIANT_MANUAL_CREATE=True, VARIANT_MANUAL_CREATE_BY_NON_ADMIN=True)
class ManualVariantEntrySampleTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="manual_variant_entry_sample_user")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def setUp(self):
        import_processing_dir = tempfile.TemporaryDirectory()
        self.addCleanup(import_processing_dir.cleanup)
        settings_override = override_settings(IMPORT_PROCESSING_DIR=import_processing_dir.name)
        settings_override.enable()
        self.addCleanup(settings_override.disable)

    def _create_manual_variants(self, sample_name=None) -> ManualVariantEntryCollection:
        with patch("annotation.manual_variant_entry.process_upload_pipeline"):
            return create_manual_variants(self.user, self.grch37, VARIANTS_TEXT, sample_name=sample_name)

    def test_default_imports_variants_only(self):
        mvec = self._create_manual_variants()
        file_upload = mvec.uploadedmanualvariantentrycollection.file_upload
        self.assertEqual(file_upload.file_type, UploadedFileTypes.MANUAL_VARIANT_ENTRY)

    def test_create_as_sample_imports_one_sample_vcf(self):
        sample_name = "my variant list"
        mvec = self._create_manual_variants(sample_name=sample_name)
        file_upload = mvec.uploadedmanualvariantentrycollection.file_upload
        self.assertEqual(file_upload.file_type, UploadedFileTypes.MANUAL_VARIANT_ENTRY_SAMPLE)

        upload_step = UploadStep.objects.create(upload_pipeline=file_upload.uploadpipeline, sort_order=0,
                                                input_filename=file_upload.path)
        ImportCreateVCFModelForManualVariantEntryTask.process_items(upload_step)

        mvec.refresh_from_db()
        vcf = mvec.vcf
        self.assertEqual(vcf.genome_build, self.grch37)  # From the upload - the VCF header has no contigs
        self.assertFalse(vcf.variant_zygosity_count)
        self.assertEqual(list(vcf.sample_set.values_list("name", flat=True)), [sample_name])
        self.assertEqual(mvec.sample.cohortsample_set.get().cohort, vcf.cohort)

    def test_create_as_sample_needs_a_name(self):
        data = {"variants_text": VARIANTS_TEXT, "genome_build": self.grch37.pk, "create_as_sample": True}
        form = ManualVariantEntryForm(data, user=self.user)
        self.assertFalse(form.is_valid())
        self.assertIn("sample_name", form.errors)

        form = ManualVariantEntryForm(data | {"sample_name": "my variant list"}, user=self.user)
        self.assertTrue(form.is_valid(), form.errors)


class ManualVariantEntryLinkingTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.user = User.objects.create_user(username="manual_variant_entry_linking_user")

    def test_each_entry_linked_across_batches(self):
        """ Two entries for the same variant each get linked, and the entry ids stay in step with
            the variants when the records span several batches """
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)

        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        coordinates = [("1", 169519049, "T", "C"), ("1", 169519049, "T", "C"), ("13", 32914438, "T", "TA")]
        expected = set()
        variant_coordinates = []
        vcf_ids = []
        for line_number, (chrom, position, ref, alt) in enumerate(coordinates, start=1):
            variant = slowly_create_test_variant(chrom, position, ref, alt, self.grch37)
            mve = ManualVariantEntry.objects.create(manual_variant_entry_collection=mvec, line_number=line_number,
                                                    entry_text=f"{chrom}:{position} {ref}>{alt}")
            expected.add((mve.pk, variant.pk))
            variant_coordinates.append(VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt))
            vcf_ids.append(mve.pk)

        vcf_filename = os.path.join(temp_dir.name, "manual_variant_entry.vcf")
        write_vcf_from_variant_coordinates(vcf_filename, variant_coordinates, vcf_ids=vcf_ids)
        file_upload = FileUpload.objects.create(path=vcf_filename, import_source=ImportSource.WEB, user=self.user,
                                                file_type=UploadedFileTypes.MANUAL_VARIANT_ENTRY)
        UploadedManualVariantEntryCollection.objects.create(file_upload=file_upload, collection=mvec)
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline, sort_order=0,
                                                input_filename=vcf_filename)

        processor = BulkManualVariantEntryLinkingVCFProcessor(upload_step, None, batch_size=2)
        for v in cyvcf2.VCF(vcf_filename):
            processor.process_entry(v)
        processor.finish()

        created = set(CreatedManualVariant.objects.filter(manual_variant_entry__manual_variant_entry_collection=mvec)
                      .values_list("manual_variant_entry_id", "variant_id"))
        self.assertEqual(created, expected)
