from datetime import timedelta

from django.dispatch import receiver

from library.health_check import HealthCheckAge, HealthCheckRequest, \
    health_check_overall_stats_signal
from sync.models import SyncRun, SyncDestination, SyncStatus


@receiver(signal=health_check_overall_stats_signal)
def sync_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    # Report when each enabled sync run was last successfully performed
    responses: list[HealthCheckAge] = []
    for sync_destination in SyncDestination.objects.all():
        last_successful_sync_run = SyncRun.objects.filter(
                destination=sync_destination,
                status__in=(SyncStatus.SUCCESS, SyncStatus.NO_RECORDS)
        ).order_by('-created').values_list('created', flat=True).first()
        responses.append(
            HealthCheckAge(
                name="SYNC (" + sync_destination.name + ")" + (" *disabled*" if not sync_destination.enabled else ""),
                now=health_request.now,
                last_performed=last_successful_sync_run,
                warning_age=timedelta(days=1)
            )
        )
    return responses
