from django.dispatch import receiver

from eventlog.models import IntegrationActivity
from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    integration_status_signal,
)


@receiver(signal=integration_status_signal)
def integration_activity_status(sender, **kwargs) -> list[IntegrationStatus]:
    """ Anything that calls IntegrationActivity.track()/record() appears on Server Status with
        no further registration - apps only write their own receiver when they want more """
    statuses: list[IntegrationStatus] = []
    for activity in IntegrationActivity.objects.order_by('name'):
        details = [
            IntegrationDetail(label="Last Run", timestamp=activity.last_attempt),
            IntegrationDetail(label="Last Success", timestamp=activity.last_success),
        ]
        if activity.last_change:
            details.append(IntegrationDetail(label="Last Change", timestamp=activity.last_change))
        if activity.last_error:
            details.append(IntegrationDetail(
                label="Last Failure",
                timestamp=activity.last_error,
                # a system failing now reads differently from one that failed once and recovered
                status="danger" if activity.failing else "warning",
                extra=activity.last_error_message
            ))
        statuses.append(IntegrationStatus(
            name=activity.name,
            details=details,
            direction=activity.direction
        ))
    return statuses
