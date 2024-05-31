import logging
import time

import celery
import requests
from django.conf import settings

from library.constants import MINUTE_SECS
from library.log_utils import report_message
from variantgrid.deployment_validation.disk_usage import get_disk_messages


@celery.shared_task
def warn_low_disk_space():
    low_disk_messages = get_disk_messages(info_messages=False)
    if low_disk_messages:
        message = "\n".join([m[1] for m in low_disk_messages])
        report_message(message=message, level='warning')


@celery.shared_task
def heartbeat():
    if settings.HEARTBEAT_URL:
        _ = requests.get(settings.HEARTBEAT_URL, timeout=MINUTE_SECS)
        # we're not overly concerned with the response


@celery.shared_task
def sleep_task(seconds: int):
    time.sleep(seconds)
    logging.info("Done sleeping...")
