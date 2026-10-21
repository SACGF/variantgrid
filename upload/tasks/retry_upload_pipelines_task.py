"""
Bulk retry of failed imports (#2188) - the upload page's "Retry failed" and the sequencing run page's
"Retry failed imports". Each retry is its own task, RETRY_STAGGER_SECONDS after the last, so relaunching a whole
run's imports doesn't recreate the memory spike that failed them.

A step whose worker was killed from outside (OOM killer, service restart) never reaches ERROR: the task_failure
handler in import_vcf_step_task only sees a child lost under a surviving master, and celery acks before running,
so the task is simply gone and the step and pipeline stay PROCESSING. Its lease (upload_step_lease) stops being
renewed though, and `fail_stuck_upload_pipelines` marks a step with a lapsed lease TIMED_OUT so a retry picks it
up - from the beat-scheduled `fail_stuck_upload_pipelines_task` and again before each retry.

Entry points: `fail_stuck_upload_pipelines_task` (beat), `fail_stuck_upload_pipelines` (before a retry),
`get_retryable_upload_pipelines` (what a user may retry), `retry_upload_pipelines` (queue them).
"""
import logging
import os

import celery
from django.contrib.auth.models import User
from django.db.models import QuerySet

from snpdb.models.models_enums import ProcessingStatus
from upload.models import UploadPipeline, UploadStep
from upload.uploaded_file_type import retry_upload_pipeline

RETRY_STAGGER_SECONDS = 30
RETRY_QUEUED_PROGRESS_STATUS = "Retry queued"
RETRYABLE_STATUSES = (ProcessingStatus.ERROR, ProcessingStatus.TIMED_OUT)
WORKER_LOST_MESSAGE = "Worker stopped renewing this step's lease (killed or restarted)"


def fail_stuck_upload_pipelines(upload_pipelines_qs: QuerySet[UploadPipeline]) -> int:
    """ Mark TIMED_OUT every running step of these pipelines whose lease has lapsed, failing the pipeline with
        it. Returns how many pipelines were failed. A step queued but not yet picked up is still CREATED and
        left alone - it is still in the broker """
    processing_qs = upload_pipelines_qs.filter(status=ProcessingStatus.PROCESSING)
    lapsed_steps_qs = UploadStep.objects.filter(UploadStep.lease_lapsed_q(), upload_pipeline__in=processing_qs) \
        .select_related("upload_pipeline")
    failed_pipeline_ids = set()
    for upload_step in lapsed_steps_qs:
        logging.warning("fail_stuck_upload_pipelines: %s lease lapsed at %s (task %s)",
                        upload_step, upload_step.lease_expires, upload_step.celery_task)
        upload_step.mark_timed_out(error_message=WORKER_LOST_MESSAGE)
        failed_pipeline_ids.add(upload_step.upload_pipeline_id)
    return len(failed_pipeline_ids)


@celery.shared_task
def fail_stuck_upload_pipelines_task():
    """ Beat-scheduled sweep so a killed import shows as TIMED_OUT without anyone clicking retry """
    if num_failed := fail_stuck_upload_pipelines(UploadPipeline.objects.all()):
        logging.warning("fail_stuck_upload_pipelines_task: timed out %d pipeline(s) whose worker died", num_failed)


def get_retryable_upload_pipelines(upload_pipelines_qs: QuerySet[UploadPipeline], user: User) -> list[UploadPipeline]:
    """ Failed pipelines the user may retry, with the uploaded file still on disk to re-read """
    retryable = []
    for upload_pipeline in upload_pipelines_qs.filter(status__in=RETRYABLE_STATUSES).select_related("file_upload__user"):
        file_upload = upload_pipeline.file_upload
        filename = file_upload.get_filename()
        if filename and os.path.exists(filename) and file_upload.can_write(user):
            retryable.append(upload_pipeline)
    return retryable


def retry_upload_pipelines(upload_pipelines: list[UploadPipeline]) -> int:
    """ Returns how many were queued. Each pipeline is claimed (failed -> CREATED) before its task is queued,
        so a second click while the first batch waits doesn't queue it again """
    num_queued = 0
    for upload_pipeline in upload_pipelines:
        claimed = UploadPipeline.objects.filter(pk=upload_pipeline.pk, status__in=RETRYABLE_STATUSES).update(
            status=ProcessingStatus.CREATED, progress_status=RETRY_QUEUED_PROGRESS_STATUS, error_summary=None)
        if claimed:
            retry_upload_pipeline_task.apply_async((upload_pipeline.pk,), countdown=num_queued * RETRY_STAGGER_SECONDS)
            num_queued += 1
    return num_queued


def retry_message(num_queued: int, num_stuck: int) -> str:
    message = f"Retrying {num_queued} failed import(s)"
    if num_stuck:
        message += f" ({num_stuck} had no worker still running them)"
    return message


@celery.shared_task
def retry_upload_pipeline_task(upload_pipeline_id: int):
    if upload_pipeline := UploadPipeline.objects.filter(pk=upload_pipeline_id).first():  # may be deleted while queued
        retry_upload_pipeline(upload_pipeline)
