from django.dispatch import receiver
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity
from snpdb.models import VCF
from variantgrid.perm_path import get_visible_url_names


@receiver(signal=health_check_signal)
def vcf_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    if get_visible_url_names().get('analyses'):
        return HealthCheckRecentActivity.simple_report(
            health_request=health_request,
            emoji=":green_book:",
            model=VCF,
            created=True)
