import os
from typing import Tuple

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import ImportSource, Sequence, md5sum_str
from upload.models import UploadedFile, UploadPipeline, UploadedVCF, UploadStep, UploadedFileTypes
from upload.vcf.bulk_genotype_vcf_processor import BulkGenotypeVCFProcessor
from upload.vcf.bulk_no_genotype_vcf_processor import BulkNoGenotypeVCFProcessor
from upload.vcf.sql_copy_files import COHORT_GENOTYPE_HEADER
from upload.vcf.vcf_import import create_vcf_from_vcf, create_cohort_genotype_collection_from_vcf


class TestVCFProcessors(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf")

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_md5_hash=md5sum_str(base), length=len(base))

    @classmethod
    def _create_fake_upload_step_and_vcf(cls, vcf_filename, vcf_reader) -> Tuple[UploadStep, UploadedVCF]:
        user = User.objects.get_or_create(username='testuser')[0]
        uploaded_file = UploadedFile.objects.create(path=vcf_filename,
                                                    import_source=ImportSource.COMMAND_LINE,
                                                    user=user,
                                                    file_type=UploadedFileTypes.VCF)

        upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                input_filename=vcf_filename,
                                                sort_order=0)
        vcf = create_vcf_from_vcf(upload_step, vcf_reader)
        get_fake_annotation_version(vcf.genome_build)
        create_cohort_genotype_collection_from_vcf(vcf, vcf_reader)
        uploaded_vcf = UploadedVCF.objects.get(upload_pipeline=upload_pipeline)
        return upload_step, uploaded_vcf

    def _test_genotype_processor(self, vcf_filename, processor_klass):
        """ I keep forgetting to adjust the columns to match the CSV """
        fast_vcf_reader = cyvcf2.VCF(vcf_filename)
        upload_step, uploaded_vcf = self._create_fake_upload_step_and_vcf(vcf_filename, fast_vcf_reader)
        cohort_genotype_collection = uploaded_vcf.vcf.cohort.cohort_genotype_collection
        processor = processor_klass(upload_step, cohort_genotype_collection, uploaded_vcf, None)

        for v in fast_vcf_reader:
            processor.process_entry(v)
            break

        cg = None
        for field_name in ["locus_cohort_genotypes", "cohort_genotypes"]:
            if f := getattr(processor, field_name, None):
                cg = f[0]
                break
        if cg is None:
            raise ValueError("Couldn't find array to retrieve cohort genotype")
        len_genotype_cols = len(cg)
        len_columns = len(COHORT_GENOTYPE_HEADER) - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS
        message = f"{processor_klass} CohortGenotypeData ({len_genotype_cols} cols) !=  CSV columns ({len_columns})"
        self.assertEqual(len_genotype_cols, len_columns, message)

    def test_no_genotype_processor(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "no_genotype.GRCh37.vcf")
        self._test_genotype_processor(vcf_filename, BulkNoGenotypeVCFProcessor)

    def test_genotype_processor(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "sample1_hg19.vcf")
        self._test_genotype_processor(vcf_filename, BulkGenotypeVCFProcessor)
