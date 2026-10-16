""" LiftoverCreateVCFTask: which of a run's AlleleLiftovers it moves to PROCESSING (#2246) """
from unittest.mock import patch

from django.test import TestCase

from library.guardian_utils import admin_bot
from snpdb.models import AlleleLiftover, LiftoverRun
from snpdb.models.models_enums import AlleleConversionTool, ImportSource, ProcessingStatus
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


class LiftoverCreateVCFTaskTest(TestCase):

    def test_precheck_errors_stay_error(self):
        grch37 = GenomeBuild.grch37()
        user = admin_bot()
        liftover = LiftoverRun.objects.create(user=user, conversion_tool=AlleleConversionTool.BCFTOOLS_LIFTOVER,
                                              genome_build=GenomeBuild.grch38(), source_genome_build=grch37)
        in_vcf = AlleleLiftover.objects.create(
            allele=create_mock_allele(slowly_create_test_variant("3", 1000, 'A', 'T', grch37), grch37),
            liftover=liftover, status=ProcessingStatus.CREATED)
        failed_precheck = AlleleLiftover.objects.create(
            allele=create_mock_allele(slowly_create_test_variant("3", 2000, 'A', 'T', grch37), grch37),
            liftover=liftover, status=ProcessingStatus.ERROR,
            error={AlleleLiftover.ERROR_JSON_MESSAGE_KEY: "Validation Error: symbolic"})

        file_upload = FileUpload.objects.create(user=user, path="/tmp/liftover.vcf", import_source=ImportSource.WEB,
                                                file_type=UploadedFileTypes.LIFTOVER)
        UploadedLiftover.objects.create(file_upload=file_upload, liftover=liftover)
        upload_step = UploadStep.objects.create(upload_pipeline=UploadPipeline.objects.create(file_upload=file_upload),
                                                name="Liftover", sort_order=1, task_type=UploadStepTaskType.CELERY)

        with patch("upload.tasks.vcf.import_vcf_tasks.bcftools_liftover", return_value=(None, 0)):
            LiftoverCreateVCFTask.process_items(upload_step)

        in_vcf.refresh_from_db()
        failed_precheck.refresh_from_db()
        self.assertEqual(in_vcf.status, ProcessingStatus.PROCESSING)
        self.assertEqual(failed_precheck.status, ProcessingStatus.ERROR)
