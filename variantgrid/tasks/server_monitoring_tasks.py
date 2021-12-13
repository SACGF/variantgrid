from typing import List, Tuple

import celery
import requests
from django.conf import settings

from library.file_utils import get_disk_usage
from library.log_utils import report_message


def get_disk_messages(directories_list: List[str] = None, info_messages=False) -> List[Tuple[str, str]]:
    if directories_list is None:
        directories_list = [settings.BASE_DIR, settings.UPLOAD_DIR, settings.ANNOTATION_VCF_DUMP_DIR]

    minimum_gigs = settings.SERVER_MIN_DISK_WARNING_GIGS
    minimum_kb = minimum_gigs * 1000000

    disk_usage = get_disk_usage()
    nice_disk_usage = get_disk_usage(human_readable=True)
    handled_mount_points = set()
    low_disk_messages = []
    for mount_point, data in disk_usage.items():
        for d in directories_list:
            if mount_point in handled_mount_points:
                continue
            if d.startswith(mount_point):
                handled_mount_points.add(mount_point)
                available = int(data["avail"])
                percent = data["percent"]
                nice_available = nice_disk_usage[mount_point]["avail"]
                message = f"Mount point '{mount_point}' ({percent} used, {nice_available} available)"
                if available > minimum_kb:
                    status = "info"
                else:
                    status = "warning"
                    message += f" is below minimum of {minimum_gigs}G"

                if status == "warning" or info_messages:
                    low_disk_messages.append((status, message))
    return low_disk_messages


@celery.shared_task
def warn_low_disk_space():
    low_disk_messages = get_disk_messages(info_messages=False)
    if low_disk_messages:
        message = "\n".join([m[1] for m in low_disk_messages])
        report_message(message=message, level='warning')


@celery.shared_task
def heartbeat():
    if settings.HEARTBEAT_URL:
        _ = requests.get(settings.HEARTBEAT_URL)
        # we're not overly concerned with the response