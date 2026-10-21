import tempfile
from datetime import timedelta
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models.models_enums import ImportSource, ProcessingStatus
from upload.enums import UploadedFileTypes
from upload.models import FileUpload, UploadPipeline, UploadStep
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
    """ A step whose worker died from outside stays PROCESSING; its lapsed lease is how the sweep notices """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="stuck_owner")

    def _make_step(self, status=ProcessingStatus.PROCESSING, lease_expires=None, started_ago=None) -> UploadStep:
        file_upload = FileUpload.objects.create(user=self.owner, name="run.vcf", path="/tmp/run.vcf",
                                                file_type=UploadedFileTypes.VCF, import_source=ImportSource.SEQAUTO)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        start_date = timezone.now() - (started_ago or timedelta())
        return UploadStep.objects.create(upload_pipeline=upload_pipeline, name="Process VCF file", script="x",
                                         status=status, start_date=start_date, lease_expires=lease_expires,
                                         celery_task="task", sort_order=0)

    def _assert_status(self, upload_step, status):
        upload_step.refresh_from_db()
        self.assertEqual(upload_step.status, status)
        self.assertEqual(upload_step.upload_pipeline.status, status)

    def test_lapsed_lease_times_out_its_pipeline(self):
        now = timezone.now()
        dead = self._make_step(lease_expires=now - timedelta(seconds=1))
        alive = self._make_step(lease_expires=now + timedelta(minutes=5))
        queued = self._make_step(status=ProcessingStatus.CREATED, lease_expires=now - timedelta(hours=1))

        self.assertEqual(fail_stuck_upload_pipelines(UploadPipeline.objects.all()), 1)

        self._assert_status(dead, ProcessingStatus.TIMED_OUT)
        self.assertEqual(dead.error_message, WORKER_LOST_MESSAGE)
        self.assertIsNotNone(dead.end_date)
        self._assert_status(alive, ProcessingStatus.PROCESSING)
        queued.refresh_from_db()
        self.assertEqual(queued.status, ProcessingStatus.CREATED)

    def test_step_without_lease_is_judged_on_start_date(self):
        """ Started before lease_expires existed """
        lease = timedelta(seconds=settings.UPLOAD_STEP_LEASE_SECONDS)
        old = self._make_step(started_ago=lease + timedelta(seconds=1))
        recent = self._make_step(started_ago=lease - timedelta(minutes=1))

        self.assertEqual(fail_stuck_upload_pipelines(UploadPipeline.objects.all()), 1)
        self._assert_status(old, ProcessingStatus.TIMED_OUT)
        self._assert_status(recent, ProcessingStatus.PROCESSING)

    def test_renew_lease_only_while_held(self):
        now = timezone.now()
        upload_step = self._make_step(lease_expires=now)

        self.assertTrue(upload_step.renew_lease("task"))
        upload_step.refresh_from_db()
        self.assertGreater(upload_step.lease_expires, now)

        self.assertFalse(upload_step.renew_lease("some other task"))
        upload_step.mark_timed_out()
        self.assertFalse(upload_step.renew_lease("task"))
