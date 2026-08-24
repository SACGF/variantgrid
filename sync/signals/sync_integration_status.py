from datetime import timedelta

from django.dispatch import receiver
from django.urls import reverse

from library.integration_status import (
    IntegrationDetail,
    IntegrationDirection,
    IntegrationStatus,
    integration_status_signal,
)
from sync.models import SyncDestination, SyncRun, SyncStatus


def _last_run_created(destination: SyncDestination, statuses=None):
    qs = SyncRun.objects.filter(destination=destination)
    if statuses:
        qs = qs.filter(status__in=statuses)
    return qs.order_by('-created').values_list('created', flat=True).first()


@receiver(signal=integration_status_signal)
def sync_integration_status(sender, **kwargs) -> list[IntegrationStatus]:
    statuses: list[IntegrationStatus] = []
    for destination in SyncDestination.objects.order_by('name'):
        last_success = _last_run_created(destination, (SyncStatus.SUCCESS, SyncStatus.NO_RECORDS))
        last_failure = _last_run_created(destination, (SyncStatus.FAILED,))

        details = [
            IntegrationDetail(label="Last Run", timestamp=_last_run_created(destination)),
            IntegrationDetail(label="Last Success", timestamp=last_success),
        ]
        if last_failure:
            details.append(IntegrationDetail(
                label="Last Failure",
                timestamp=last_failure,
                status="danger" if not last_success or last_failure > last_success else "warning"
            ))

        statuses.append(IntegrationStatus(
            name=f"Sync ({destination.name})",
            details=details,
            direction=IntegrationDirection.OUTBOUND if destination.config.get("direction") == "upload" else IntegrationDirection.INBOUND,
            url=reverse('admin:sync_syncdestination_change', args=[destination.pk]),
            record_count=SyncRun.objects.filter(destination=destination).count(),
            enabled=destination.enabled,
            warning_age=timedelta(days=1)
        ))
    return statuses
