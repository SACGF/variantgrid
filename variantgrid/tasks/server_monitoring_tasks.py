from django.conf import settings

from library.file_utils import get_disk_usage
from library.log_utils import report_message
from variantgrid import celery
import requests


@celery.task
def warn_low_disk_space():
    DIRECTORIES_WE_CARE_ABOUT = [settings.BASE_DIR]
    minimum_gigs = settings.ROLLBAR_MIN_DISK_WARNING_GIGS
    minimum_kb = minimum_gigs * 1000000

    disk_usage = get_disk_usage()
    nice_disk_usage = get_disk_usage(human_readable=True)
    low_disk_messages = []
    for mount_point, data in disk_usage.items():
        for d in DIRECTORIES_WE_CARE_ABOUT:
            if d.startswith(mount_point):
                available = int(data["avail"])
                if available < minimum_kb:
                    percent = data["percent"]
                    nice_available = nice_disk_usage[mount_point]["avail"]
                    params = (mount_point, percent, nice_available, minimum_gigs)
                    message = "Mount point '%s' (%s used, %s available) is below minimum of %sG" % params
                    low_disk_messages.append(message)

    message = "\n".join(low_disk_messages)
    report_message(message=message, level='warning')


@celery.task
def heartbeat():
    if settings.HEARTBEAT_URL:
        r = requests.get(settings.HEARTBEAT_URL)
        # we're not overly concerned with the response