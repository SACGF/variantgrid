from django.db import models
from django.utils import timezone


class BeaconInboundQuery(models.Model):
    """ Audit row for one inbound Beacon /g_variants query: the exact request we were
        sent (so we can see what people search for) plus a small response summary for
        health-check aggregation. Mirrors mme.models.MMEInboundQuery. """
    created = models.DateTimeField(default=timezone.now)
    request_json = models.JSONField()                   # exact Beacon query we received
    granularity = models.CharField(max_length=8)        # boolean | count | record (returned)
    authenticated = models.BooleanField(default=False)  # anon (public tier) vs authed requester
    observations_exists = models.BooleanField(default=False)
    observations_count = models.IntegerField(default=0)
    classifications_exists = models.BooleanField(default=False)
    classifications_count = models.IntegerField(default=0)

    def __str__(self):
        return f"BeaconInboundQuery({self.created}, granularity={self.granularity})"


class BeaconQueryCache(models.Model):
    """ Cached outbound result for one (variant, node): presence/count read back from an
        external Beacon, refreshed live on expiry (settings.BEACON_QUERY_CACHE_DAYS). §9.5.
        Only a public coordinate is ever sent to reach these results. """
    variant = models.ForeignKey("snpdb.Variant", on_delete=models.CASCADE)
    node_id = models.CharField(max_length=64)           # key into settings.BEACON_QUERY_NODES
    exists = models.BooleanField(null=True, blank=True)  # None when the node didn't answer
    count = models.IntegerField(null=True, blank=True)  # None when the node returns no count
    error = models.TextField(null=True, blank=True)     # node unreachable / errored -> "unavailable"
    created = models.DateTimeField(default=timezone.now)

    class Meta:
        unique_together = ("variant", "node_id")

    def __str__(self):
        return f"BeaconQueryCache(variant={self.variant_id}, node={self.node_id}, exists={self.exists})"
