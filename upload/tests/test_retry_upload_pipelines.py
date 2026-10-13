import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models.models_enums import ImportSource, ProcessingStatus
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline
from upload.tasks.retry_upload_pipelines_task import (
    RETRY_STAGGER_SECONDS,
    get_retryable_upload_pipelines,
    retry_upload_pipelines,
)


class RetryUploadPipelinesTest(TestCase):
    """ Bulk retry of failed imports (#2188) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="bulk_retry_owner")
        cls.other_user = User.objects.create_user(username="bulk_retry_other")

    def setUp(self):
        self.vcf_file = tempfile.NamedTemporaryFile(suffix=".vcf")

    def tearDown(self):
        self.vcf_file.close()

    def _make_pipeline(self, status=ProcessingStatus.ERROR, path=None) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.owner, name="run.vcf", path=path or self.vcf_file.name,
                                                file_type=UploadedFileTypes.VCF, import_source=ImportSource.SEQAUTO)
        return UploadPipeline.objects.create(status=status, file_upload=file_upload)

    def test_retryable_only_failed_writable_and_on_disk(self):
        failed = self._make_pipeline()
        self._make_pipeline(status=ProcessingStatus.SUCCESS)
        self._make_pipeline(path="/no/such/file.vcf")
        qs = UploadPipeline.objects.all()

        self.assertEqual(get_retryable_upload_pipelines(qs, self.owner), [failed])
        self.assertEqual(get_retryable_upload_pipelines(qs, self.other_user), [])

    @patch("upload.tasks.retry_upload_pipelines_task.retry_upload_pipeline_task.apply_async")
    def test_retries_are_staggered_and_queued_once(self, apply_async):
        upload_pipelines = [self._make_pipeline(), self._make_pipeline()]

        self.assertEqual(retry_upload_pipelines(upload_pipelines), 2)
        countdowns = [call.kwargs["countdown"] for call in apply_async.call_args_list]
        self.assertEqual(countdowns, [0, RETRY_STAGGER_SECONDS])
        self.assertEqual(retry_upload_pipelines(upload_pipelines), 0)  # second click while queued
