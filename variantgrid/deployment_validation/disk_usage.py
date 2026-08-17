""" Deployment view of disk space: which directories this site cares about, and the minimum it wants
    free on them. The measurement itself is library.utils.file_utils.DiskUsage. """
from typing import Optional

from django.conf import settings

from library.utils.file_utils import DiskUsage, get_disk_usage


def get_disk_usage_objects(directories_list: Optional[list[str]] = None) -> list[DiskUsage]:
    """ Mount points serving the given directories, for operator-facing reporting.

        Deliberately broad: every mount point that is a textual prefix of a directory is reported, so a
        directory on '/data' yields both '/data' and '/'. That is wanted here - '/' holds the database and
        the site, so it must stay watched no matter which deeper filesystem a directory sits on.

        A caller that needs the one filesystem a path actually lives on (eg the annotation dispatcher
        asking whether there is room for another VEP dump) wants get_disk_usage_for_directory, whose
        longest match resolves '/data/annotation_scratch' to '/data' alone. """
    if directories_list is None:
        directories_list = [settings.BASE_DIR, settings.UPLOAD_DIR, settings.ANNOTATION_VCF_DUMP_DIR]

    disk_usage = get_disk_usage()
    nice_disk_usage = get_disk_usage(human_readable=True)
    handled_mount_points = set()
    disk_usages: list[DiskUsage] = []
    for mount_point, data in disk_usage.items():
        for d in directories_list:
            if mount_point in handled_mount_points:
                continue
            if d.startswith(mount_point):
                handled_mount_points.add(mount_point)
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
