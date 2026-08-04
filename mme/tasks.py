from variantgrid.celery import app

from mme.client import submit
from mme.metrics import get_metrics
from mme.models import MMESubmission
from mme.notifications import notify_depositors


@app.task(queue='web_workers')
def submit_mme_submission_task(submission_id: int) -> None:
    """ Build the profile, POST to the remote node and persist results.
        Runs on web_workers (network-bound) so the UI POST returns immediately. """
    submission = MMESubmission.objects.get(pk=submission_id)
    submit(submission)


@app.task(queue='web_workers')
def notify_mme_depositors_task(inbound_query_id: int) -> None:
    """ Depositor notification, off the inbound request's critical path. """
    notify_depositors(inbound_query_id)


@app.task(queue='db_workers')
def refresh_mme_metrics_task() -> None:
    """ Recompute and cache the MME Metrics API payload. """
    get_metrics(refresh=True)
