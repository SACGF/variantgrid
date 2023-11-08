from django.dispatch import receiver

from library.health_check import HealthCheckRequest, HealthCheckCapacity, \
    health_check_overall_stats_signal
from variantgrid.tasks.server_monitoring_tasks import get_disk_usage_objects


@receiver(signal=health_check_overall_stats_signal)
def disk_usage_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    checks = []
    for disk_usage in get_disk_usage_objects():
        checks.append(HealthCheckCapacity(
            name=f"Mount Point '{disk_usage.mount_point}",
            used=disk_usage.percent_nice,
            available=disk_usage.available_nice,
            warning=not disk_usage.has_safe_capacity
        ))
    if not checks:
        checks.append(HealthCheckCapacity(
            name="Disk Usage Unknown",
            used="?",
            available="?"
        ))
    return checks
