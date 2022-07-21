from django.dispatch import receiver

from analysis.models import Analysis
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity
from variantgrid.perm_path import get_visible_url_names


@receiver(signal=health_check_signal)
def analysis_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    if get_visible_url_names().get('analyses'):
        return HealthCheckRecentActivity.simple_report(
            health_request=health_request,
            emoji=":orange_book:",
            model=Analysis,
            created=True,
            modified=True)
