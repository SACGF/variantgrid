""" Deployment view of disk space: which directories this site cares about, and the minimum it wants
    free on them. The measurement itself is library.utils.file_utils.DiskUsage. """
from typing import Optional

from django.conf import settings

from library.utils.file_utils import DiskUsage, get_disk_usage, get_mount_point_for_directory


def get_disk_usage_objects(directories_list: Optional[list[str]] = None) -> list[DiskUsage]:
    """ One DiskUsage per distinct mount point serving the given directories. """
    if directories_list is None:
        directories_list = [settings.BASE_DIR, settings.UPLOAD_DIR, settings.ANNOTATION_VCF_DUMP_DIR]

    disk_usage = get_disk_usage()
    nice_disk_usage = get_disk_usage(human_readable=True)
    handled_mount_points = set()
    disk_usages: list[DiskUsage] = []
    for d in directories_list:
        mount_point = get_mount_point_for_directory(d, disk_usage)
        if mount_point is None or mount_point in handled_mount_points:
            continue
        handled_mount_points.add(mount_point)
        data = disk_usage[mount_point]
        disk_usages.append(
            DiskUsage(
                mount_point=mount_point,
                available_kb=int(data["avail"]),
                percent_nice=data["percent"],
                available_nice=nice_disk_usage[mount_point]["avail"]
            )
        )
    return disk_usages


def get_disk_messages(directories_list: list[str] = None, info_messages=False,
                      min_gigs: Optional[float] = None) -> list[tuple[str, str]]:
    if min_gigs is None:
        min_gigs = settings.SERVER_MIN_DISK_WARNING_GIGS

    disk_usages = get_disk_usage_objects(directories_list)
    if not info_messages:
        disk_usages = [du for du in disk_usages if not du.has_capacity(min_gigs)]

    return [du.as_status_message(min_gigs) for du in disk_usages]
