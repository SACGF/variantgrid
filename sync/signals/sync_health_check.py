from datetime import timedelta
from typing import List

from django.dispatch import receiver

from library.health_check import health_check_signal, HealthCheckAge, HealthCheckRequest
from sync.models import SyncRun, SyncDestination, SyncStatus


@receiver(signal=health_check_signal)
def sync_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    # Report when each enabled sync run was last successfully performed
    responses: List[HealthCheckAge] = []
    for sync_destination in SyncDestination.objects.filter(enabled=True):
        last_successful_sync_run = SyncRun.objects.filter(
                destination=sync_destination,
                status__in=(SyncStatus.SUCCESS, SyncStatus.NO_RECORDS)
        ).order_by('-created').values_list('created', flat=True).first()
        responses.append(
            HealthCheckAge(
                name=sync_destination.name,
                now=health_request.now,
                last_performed=last_successful_sync_run,
                warning_age=timedelta(days=1)
            )
        )
    return responses
