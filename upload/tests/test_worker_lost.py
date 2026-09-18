"""
A VCF import step whose worker child is killed (SIGTERM/SIGKILL, or revoke with terminate) never runs its own
error handling, so the celery master's task_failure / task_revoked receivers have to fail it instead.
"""
from types import SimpleNamespace

from celery.exceptions import WorkerLostError
from celery.signals import task_failure, task_revoked
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models.models_enums import ImportSource
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
    UploadStepTaskType,
)
from upload.tasks.vcf.import_vcf_tasks import PreprocessVCFTask


class WorkerLostTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="worker_lost_user", password="x")

    def _make_step(self, **kwargs) -> UploadStep:
        file_upload = FileUpload.objects.create(user=self.user, name="tso500.vcf", path="/tmp/tso500.vcf",
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.WEB)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING,
                                                        file_upload=file_upload)
        return UploadStep.objects.create(upload_pipeline=upload_pipeline, name=UploadStep.PREPROCESS_VCF_NAME,
                                         sort_order=1, task_type=UploadStepTaskType.CELERY,
                                         status=ProcessingStatus.PROCESSING, start_date=timezone.now(), **kwargs)

    def _assert_failed(self, upload_step: UploadStep, error_text: str):
        upload_step.refresh_from_db()
        self.assertEqual(upload_step.status, ProcessingStatus.ERROR)
        self.assertIn(error_text, upload_step.error_message)
        self.assertIsNotNone(upload_step.end_date)
        self.assertEqual(upload_step.upload_pipeline.status, ProcessingStatus.ERROR)

    def test_worker_lost_fails_step_and_pipeline(self):
        upload_step = self._make_step()
        task_failure.send(sender=PreprocessVCFTask, task_id="x", args=(upload_step.pk, 0), kwargs={},
                          exception=WorkerLostError("Worker exited prematurely: signal 15 (SIGTERM)"))
        self._assert_failed(upload_step, "signal 15")

    def test_terminated_fails_step_and_pipeline(self):
        upload_step = self._make_step()
        task_revoked.send(sender=PreprocessVCFTask, request=SimpleNamespace(args=(upload_step.pk, 0)),
                          terminated=True, signum=15, expired=False)
        self._assert_failed(upload_step, "terminated")

    def test_finished_step_left_alone(self):
        """ The child can die after run() has already closed the step """
        upload_step = self._make_step(end_date=timezone.now())
        UploadStep.objects.filter(pk=upload_step.pk).update(status=ProcessingStatus.SUCCESS)
        task_failure.send(sender=PreprocessVCFTask, task_id="x", args=(upload_step.pk, 0), kwargs={},
                          exception=WorkerLostError("Worker exited prematurely: signal 15 (SIGTERM)"))
        upload_step.refresh_from_db()
        self.assertEqual(upload_step.status, ProcessingStatus.SUCCESS)
        self.assertEqual(upload_step.upload_pipeline.status, ProcessingStatus.PROCESSING)
