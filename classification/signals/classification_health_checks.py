from typing import List

from django.dispatch import receiver

from classification.models import Classification
from library.health_check import health_check_signal, \
    HealthCheckRequest, HealthCheckTotalAmount, HealthCheckRecentActivity


@receiver(signal=health_check_signal)
def classifications_health_check_count(sender, health_request, **kwargs):
    total_shared = Classification.dashboard_total_shared_classifications()
    total_unshared = Classification.dashboard_total_unshared_classifications()
    total = total_unshared + total_shared
    if total:
        percent_shared = 100.0 * float(total_shared) / float(total)

        return HealthCheckTotalAmount(
            emoji=":blue_book:",
            amount=total,
            name="Classifications",
            extra=f"{int(percent_shared)}% shared"
        )


@receiver(signal=health_check_signal)
def classification_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    classifications_of_interest = Classification.dashboard_report_classifications_of_interest(since=health_request.since)
    new_classification_count = Classification.dashboard_report_new_classifications(since=health_request.since)

    return [
        HealthCheckRecentActivity(
            emoji=":blue_book:",
            amount=new_classification_count,
            name="Classifications",
            sub_type="Created"
        ),
        HealthCheckRecentActivity(
            emoji=":blue_book:",
            amount=len(classifications_of_interest),
            name="Classifications",
            sub_type="Of Interest"
        )
    ]