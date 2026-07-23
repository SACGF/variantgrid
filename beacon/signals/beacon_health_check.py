from django.conf import settings
from django.db.models import Count
from django.dispatch import receiver

from beacon.models import BeaconInboundQuery
from library.health_check import HealthCheckRequest, HealthCheckRecentActivity, \
    health_check_overall_stats_signal


@receiver(signal=health_check_overall_stats_signal)
def beacon_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    """ Nightly Slack line: how many inbound Beacon queries arrived in the window, split by
        granularity, plus how many hit each dataset. Gated on BEACON_ENABLED so it only
        appears where Beacon is actually on. """
    if not settings.BEACON_ENABLED:
        return []

    qs = BeaconInboundQuery.objects.filter(created__gte=health_request.since,
                                           created__lt=health_request.now)
    total = qs.count()
    by_gran = dict(qs.values_list("granularity").annotate(n=Count("pk")))
    obs_hits = qs.filter(observations_exists=True).count()
    cls_hits = qs.filter(classifications_exists=True).count()
    return [HealthCheckRecentActivity(
        emoji="📡", name="Beacon queries", amount=total, sub_type="received",
        extra=(f"granularity {by_gran}; hits — observations {obs_hits}, "
               f"classifications {cls_hits}"),
        stand_alone=True,
    )]
