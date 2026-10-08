from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase

from snpdb.models.models_enums import ImportSource
from snpdb.views.datatable_view import datatable_response
from upload.grids import FileUploadColumns
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
    UploadStepTaskType,
)
from upload.vcf.vcf_ref_check import VCFRefMismatchError


class UploadUserErrorTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="upload_user_error_user")

    def _fail_step(self, exception: Exception) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.user, name="wrong_build.vcf", path="/tmp/wrong_build.vcf",
                                                file_type=UploadedFileTypes.VCF, import_source=ImportSource.WEB)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline, name=UploadStep.PREPROCESS_VCF_NAME,
                                                sort_order=1, task_type=UploadStepTaskType.CELERY,
                                                status=ProcessingStatus.PROCESSING)
        with patch("upload.models.models.report_exc_info") as report_exc_info:
            try:
                raise exception
            except Exception as e:
                upload_step.error_exception(e)
        self.reported = report_exc_info.called
        upload_pipeline.refresh_from_db()
        return upload_pipeline

    def test_user_error_has_summary_and_no_traceback(self):
        error = VCFRefMismatchError("SNV REF bases: 75% mismatch", summary="Wrong genome build - looks like GRCh38")
        upload_pipeline = self._fail_step(error)
        self.assertEqual(upload_pipeline.status, ProcessingStatus.ERROR)
        self.assertEqual(upload_pipeline.error_summary, "Wrong genome build - looks like GRCh38")
        self.assertEqual(upload_pipeline.progress_status, "Error: SNV REF bases: 75% mismatch")
        self.assertFalse(self.reported)

        request = RequestFactory().get("/")
        request.user = self.user
        status = datatable_response(FileUploadColumns(request))["data"][0]["status"]
        self.assertEqual(status["summary"], "Wrong genome build - looks like GRCh38")
        self.assertEqual(status["title"], "Error: SNV REF bases: 75% mismatch")

    def test_other_error_keeps_traceback(self):
        upload_pipeline = self._fail_step(ValueError("bug"))
        self.assertIsNone(upload_pipeline.error_summary)
        self.assertIn("Traceback", upload_pipeline.progress_status)
        self.assertTrue(self.reported)
