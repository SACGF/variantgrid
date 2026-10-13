"""
Lease heartbeat for a running UploadStep (#2188). The worker renews `UploadStep.lease_expires` on a background
thread for as long as the step's task runs, so a step whose lease has lapsed had its worker killed from outside
(OOM killer, service restart) - the one failure celery's own signals never report, as it acks before running.
`fail_stuck_upload_pipelines` (retry_upload_pipelines_task) marks those TIMED_OUT.

Entry point: `UploadStepLeaseHeartbeat`, a context manager around `ImportVCFStepTask.process_items`.
"""
import logging
import threading

from django.conf import settings
from django.db import connection

from library.log_utils import log_traceback
from upload.models.models import UploadStep


class UploadStepLeaseHeartbeat:
    """ Renews the step's lease every UPLOAD_STEP_LEASE_HEARTBEAT_SECONDS while the task runs. The renew is
        guarded on our own celery task id: once the step has been marked timed out and retried under another
        task the update matches nothing and the heartbeat stops, so this worker never freshens a lease the
        retry now owns """

    def __init__(self, upload_step: UploadStep, celery_task: str):
        self.upload_step = upload_step
        self.celery_task = celery_task
        self.interval = settings.UPLOAD_STEP_LEASE_HEARTBEAT_SECONDS
        self._stop = threading.Event()
        self._thread = None

    def _run(self):
        try:
            while not self._stop.wait(self.interval):
                try:
                    if not self.upload_step.renew_lease(self.celery_task):
                        logging.warning("Lease heartbeat: %s no longer held by task %s - stopping",
                                        self.upload_step, self.celery_task)
                        self._stop.set()
                except Exception:
                    log_traceback()  # a DB blip costs one tick, not the step
        finally:
            connection.close()

    def __enter__(self):
        self._thread = threading.Thread(target=self._run, daemon=True,
                                        name=f"upload-step-lease-{self.upload_step.pk}")
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join(timeout=self.interval)
        return False
