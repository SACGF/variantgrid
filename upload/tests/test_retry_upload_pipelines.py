import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models.models_enums import ImportSource, ProcessingStatus
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline, UploadStep
from upload.tasks.retry_upload_pipelines_task import (
    RETRY_STAGGER_SECONDS,
    WORKER_LOST_MESSAGE,
    fail_stuck_upload_pipelines,
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

    def test_timed_out_is_retryable(self):
        timed_out = self._make_pipeline(status=ProcessingStatus.TIMED_OUT)
        self.assertEqual(get_retryable_upload_pipelines(UploadPipeline.objects.all(), self.owner), [timed_out])


class FailStuckUploadPipelinesTest(TestCase):
    """ A step whose worker died from outside stays PROCESSING until the retry notices no worker holds its task """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="stuck_owner")

    def _make_processing_step(self, celery_task) -> UploadStep:
        file_upload = FileUpload.objects.create(user=self.owner, name="run.vcf", path="/tmp/run.vcf",
                                                file_type=UploadedFileTypes.VCF, import_source=ImportSource.SEQAUTO)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        return UploadStep.objects.create(upload_pipeline=upload_pipeline, name="Process VCF file", script="x",
                                         status=ProcessingStatus.PROCESSING, start_date=timezone.now(),
                                         celery_task=celery_task, sort_order=0)

    @patch("upload.tasks.retry_upload_pipelines_task._live_celery_task_ids", return_value={"alive"})
    def test_step_with_no_worker_times_out_its_pipeline(self, _live):
        dead = self._make_processing_step("dead")
        alive = self._make_processing_step("alive")
        queued = self._make_processing_step(None)  # still in the broker

        self.assertEqual(fail_stuck_upload_pipelines(UploadPipeline.objects.all()), 1)

        dead.refresh_from_db()
        self.assertEqual(dead.status, ProcessingStatus.TIMED_OUT)
        self.assertEqual(dead.error_message, WORKER_LOST_MESSAGE)
        self.assertIsNotNone(dead.end_date)
        self.assertEqual(dead.upload_pipeline.status, ProcessingStatus.TIMED_OUT)
        for step in (alive, queued):
            step.refresh_from_db()
            self.assertEqual(step.status, ProcessingStatus.PROCESSING)
            self.assertEqual(step.upload_pipeline.status, ProcessingStatus.PROCESSING)

    @patch("upload.tasks.retry_upload_pipelines_task._live_celery_task_ids", return_value=None)
    def test_nothing_fails_when_no_worker_answers(self, _live):
        step = self._make_processing_step("unknown")
        self.assertEqual(fail_stuck_upload_pipelines(UploadPipeline.objects.all()), 0)
        step.refresh_from_db()
        self.assertEqual(step.status, ProcessingStatus.PROCESSING)
