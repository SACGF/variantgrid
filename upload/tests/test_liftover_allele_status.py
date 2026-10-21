""" A liftover run's AlleleLiftovers reach a finished status however its upload pipeline ends (#2246) """
from unittest.mock import patch

from django.test import TestCase

from library.guardian_utils import admin_bot
from snpdb.enums import AlleleConversionTool, ImportSource, ProcessingStatus
from snpdb.models import AlleleLiftover, LiftoverRun
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant
from upload.models import (
    FileUpload,
    UploadedFileTypes,
    UploadedLiftover,
    UploadPipeline,
    UploadStep,
    UploadStepTaskType,
)
from upload.tasks.vcf.import_vcf_tasks import LiftoverCreateVCFTask


class LiftoverAlleleStatusTest(TestCase):

    def setUp(self):
        grch37 = GenomeBuild.grch37()
        user = admin_bot()
        liftover = LiftoverRun.objects.create(user=user, conversion_tool=AlleleConversionTool.BCFTOOLS_LIFTOVER,
                                              genome_build=GenomeBuild.grch38(), source_genome_build=grch37)

        def _allele_liftover(position, status, error=None):
            allele = create_mock_allele(slowly_create_test_variant("3", position, 'A', 'T', grch37), grch37)
            return AlleleLiftover.objects.create(allele=allele, liftover=liftover, status=status, error=error)

        self.in_vcf = _allele_liftover(1000, ProcessingStatus.CREATED)
        self.failed_precheck = _allele_liftover(2000, ProcessingStatus.ERROR,
                                                {AlleleLiftover.ERROR_JSON_MESSAGE_KEY: "Validation Error: symbolic"})
        file_upload = FileUpload.objects.create(user=user, path="/tmp/liftover.vcf", import_source=ImportSource.WEB,
                                                file_type=UploadedFileTypes.LIFTOVER)
        UploadedLiftover.objects.create(file_upload=file_upload, liftover=liftover)
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload, status=ProcessingStatus.PROCESSING)
        self.upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline, name="Liftover", sort_order=1,
                                                     task_type=UploadStepTaskType.CELERY)

    def _assert_statuses(self, in_vcf_status):
        self.in_vcf.refresh_from_db()
        self.failed_precheck.refresh_from_db()
        self.assertEqual(self.in_vcf.status, in_vcf_status)
        self.assertEqual(self.failed_precheck.status, ProcessingStatus.ERROR)
        self.assertEqual(self.failed_precheck.error[AlleleLiftover.ERROR_JSON_MESSAGE_KEY], "Validation Error: symbolic")

    def test_create_vcf_leaves_precheck_errors(self):
        with patch("upload.tasks.vcf.import_vcf_tasks.bcftools_liftover", return_value=(None, 0)):
            LiftoverCreateVCFTask.process_items(self.upload_step)
        self._assert_statuses(ProcessingStatus.PROCESSING)

    def test_timed_out_pipeline_errors_running_alleles(self):
        AlleleLiftover.objects.filter(pk=self.in_vcf.pk).update(status=ProcessingStatus.PROCESSING)
        self.upload_step.mark_timed_out()
        self._assert_statuses(ProcessingStatus.ERROR)
