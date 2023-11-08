from django.db.models import Max
from django.dispatch import receiver

from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckRecentActivity, \
    health_check_overall_stats_signal, HealthCheckTotalAmount
from snpdb.models import VCF, Variant, Sample
from variantgrid.perm_path import get_visible_url_names


@receiver(signal=health_check_signal)
def vcf_health_check_activity(sender, health_request: HealthCheckRequest, **kwargs):
    if get_visible_url_names().get('analyses'):
        return HealthCheckRecentActivity.simple_report(
            health_request=health_request,
            emoji=":green_book:",
            model=VCF,
            created=True)


@receiver(signal=health_check_overall_stats_signal)
def vcf_health_check_overall_stats(sender, **kwargs):
    return HealthCheckTotalAmount(
        emoji=":green_book:",
        amount=VCF.objects.count(),
        name="VCFs"
    )


@receiver(signal=health_check_overall_stats_signal)
def max_variants_overall_check(sender, **kwargs):
    max_variants = Variant.objects.all().aggregate(max_pk=Max("pk"))
    return HealthCheckTotalAmount(
        emoji=":green_book:",
        amount=max_variants["max_pk"],
        name="Max Variants"
    )


@receiver(signal=health_check_overall_stats_signal)
def sample_overall_check(sender, **kwargs):
    return HealthCheckTotalAmount(
        emoji=":green_book:",
        amount=Sample.objects.count(),
        name="Samples"
    )