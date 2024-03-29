import celery

from library.log_utils import report_message


@celery.shared_task(track_started=True)
def validate_alleles():
    report_message("Running automated validation - currently nothing to do", level="info")
    pass
