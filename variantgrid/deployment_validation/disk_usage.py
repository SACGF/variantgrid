from dataclasses import dataclass
from django.conf import settings
from typing import Optional
from library.utils.file_utils import get_disk_usage


@dataclass
class DiskUsage:
    mount_point: str
    available_kb: int
    # human readable below
    available_nice: str
    percent_nice: str

    @property
    def has_safe_capacity(self):
        return self.available_kb >= settings.SERVER_MIN_DISK_WARNING_GIGS * 1000000

    @property
    def as_status_message(self) -> tuple[str, str]:
        message = f"Mount point '{self.mount_point}' ({self.percent_nice} used, {self.available_nice} available)"
        if self.has_safe_capacity:
            status = "info"
        else:
            status = "warning"
            message += f" is below minimum of {settings.SERVER_MIN_DISK_WARNING_GIGS}G"

        return status, message


def get_disk_usage_objects(directories_list: Optional[list[str]] = None) -> list[DiskUsage]:
    if directories_list is None:
        directories_list = [settings.BASE_DIR, settings.UPLOAD_DIR, settings.ANNOTATION_VCF_DUMP_DIR]

    minimum_gigs = settings.SERVER_MIN_DISK_WARNING_GIGS
    minimum_kb = minimum_gigs * 1000000

    disk_usage = get_disk_usage()
    nice_disk_usage = get_disk_usage(human_readable=True)
    handled_mount_points = set()
    low_disk_messages = []
    disk_usages: list[DiskUsage] = []
    for mount_point, data in disk_usage.items():
        for d in directories_list:
            if mount_point in handled_mount_points:
                continue
            if d.startswith(mount_point):
                handled_mount_points.add(mount_point)
                available = int(data["avail"])
                percent = data["percent"]
                nice_available = nice_disk_usage[mount_point]["avail"]

                disk_usages.append(
                    DiskUsage(
                        mount_point=mount_point,
                        available_kb=available,
                        percent_nice=percent,
                        available_nice=nice_available
                    )
                )
    return disk_usages


def get_disk_messages(directories_list: list[str] = None, info_messages=False) -> list[tuple[str, str]]:
    disk_usages = get_disk_usage_objects(directories_list)
    if not info_messages:
        disk_usages = [du for du in disk_usages if not du.has_safe_capacity]

    return [du.as_status_message for du in disk_usages]
