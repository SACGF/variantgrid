from django.dispatch import receiver

from library.health_check import HealthCheckAge, HealthCheckRequest, health_check_overall_stats_signal
from library.integration_status import get_integration_statuses


@receiver(signal=health_check_overall_stats_signal)
def integration_health_check(sender, health_request: HealthCheckRequest, **kwargs) -> list[HealthCheckAge]:
    """ An integration that declared a warning_age reaches the nightly digest when it goes quiet """
    statuses, _ = get_integration_statuses()
    return [
        HealthCheckAge(
            name=status.name,
            now=health_request.now,
            last_performed=status.last_activity,
            warning_age=status.warning_age
        )
        for status in statuses if status.enabled and status.warning_age
    ]
