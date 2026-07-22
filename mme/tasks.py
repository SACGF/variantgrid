from variantgrid.celery import app

from mme.client import submit
from mme.models import MMESubmission


@app.task(queue='web_workers')
def submit_mme_submission_task(submission_id: int) -> None:
    """ Build the profile, POST to the remote node and persist results.
        Runs on web_workers (network-bound) so the UI POST returns immediately. """
    submission = MMESubmission.objects.get(pk=submission_id)
    submit(submission)
