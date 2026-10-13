"""
Bulk retry of failed imports (#2188) - the upload page's "Retry failed" and the sequencing run page's
"Retry failed imports". Each retry is its own task, RETRY_STAGGER_SECONDS after the last, so relaunching a whole
run's imports doesn't recreate the memory spike that failed them.

Entry points: `get_retryable_upload_pipelines` (what a user may retry), `retry_upload_pipelines` (queue them).
"""
import os

import celery
from django.contrib.auth.models import User
from django.db.models import QuerySet

from snpdb.models.models_enums import ProcessingStatus
from upload.models import UploadPipeline
from upload.uploaded_file_type import retry_upload_pipeline

RETRY_STAGGER_SECONDS = 30
RETRY_QUEUED_PROGRESS_STATUS = "Retry queued"


def get_retryable_upload_pipelines(upload_pipelines_qs: QuerySet[UploadPipeline], user: User) -> list[UploadPipeline]:
    """ Failed pipelines the user may retry, with the uploaded file still on disk to re-read """
    retryable = []
    for upload_pipeline in upload_pipelines_qs.filter(status=ProcessingStatus.ERROR).select_related("file_upload__user"):
        file_upload = upload_pipeline.file_upload
        filename = file_upload.get_filename()
        if filename and os.path.exists(filename) and file_upload.can_write(user):
            retryable.append(upload_pipeline)
    return retryable


def retry_upload_pipelines(upload_pipelines: list[UploadPipeline]) -> int:
    """ Returns how many were queued. Each pipeline is claimed (ERROR -> CREATED) before its task is queued,
        so a second click while the first batch waits doesn't queue it again """
    num_queued = 0
    for upload_pipeline in upload_pipelines:
        claimed = UploadPipeline.objects.filter(pk=upload_pipeline.pk, status=ProcessingStatus.ERROR).update(
            status=ProcessingStatus.CREATED, progress_status=RETRY_QUEUED_PROGRESS_STATUS, error_summary=None)
        if claimed:
            retry_upload_pipeline_task.apply_async((upload_pipeline.pk,), countdown=num_queued * RETRY_STAGGER_SECONDS)
            num_queued += 1
    return num_queued


@celery.shared_task
def retry_upload_pipeline_task(upload_pipeline_id: int):
    if upload_pipeline := UploadPipeline.objects.filter(pk=upload_pipeline_id).first():  # may be deleted while queued
        retry_upload_pipeline(upload_pipeline)
