""" Variant-page 'External Beacons' section (§9.4).

A normal HTML-fragment view (not the DRF inbound endpoint) lazy-loaded into the variant
details page. Reads the (variant, node) cache; for stale/missing nodes it fans out the
outbound query live, upserts the cache, and renders a small node -> exists/count table.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import timedelta

from django.conf import settings
from django.http import HttpResponse
from django.shortcuts import render
from django.utils import timezone

from beacon.client import query_node
from beacon.models import BeaconQueryCache
from beacon.query_targets import evaluate_queries
from snpdb.models import Variant, GenomeBuild


def external_beacons(request, variant_id: int, genome_build_name: str):
    """ Render the external-Beacon presence/count table for a variant. Gated by
        BEACON_OUTBOUND_ENABLED, then per node: a variant only fans out to servers whose
        domain it matches (§9). Nodes that don't apply are still listed with the reason they
        were skipped, so the page shows what could load and why the rest didn't. """
    if not settings.BEACON_OUTBOUND_ENABLED:
        return HttpResponse("")

    variant = Variant.objects.get(pk=variant_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    evaluations = evaluate_queries(variant, genome_build, settings.BEACON_QUERY_NODES)
    eligible = {e.node_id: e.params for e in evaluations if e.eligible}
    skipped = [e for e in evaluations if not e.eligible]

    fresh = {}
    if eligible:
        node_ids = list(eligible)
        cutoff = timezone.now() - timedelta(days=settings.BEACON_QUERY_CACHE_DAYS)
        fresh = {c.node_id: c for c in BeaconQueryCache.objects.filter(
            variant=variant, node_id__in=node_ids, created__gte=cutoff)}

        stale_nodes = [nid for nid in node_ids if nid not in fresh]
        if stale_nodes:
            with ThreadPoolExecutor(max_workers=len(stale_nodes)) as executor:
                results = list(executor.map(lambda nid: query_node(nid, eligible[nid]), stale_nodes))
            for result in results:
                cache, _ = BeaconQueryCache.objects.update_or_create(
                    variant=variant, node_id=result["node_id"],
                    defaults={"exists": result["exists"], "count": result["count"],
                              "error": result["error"], "created": timezone.now()})
                fresh[result["node_id"]] = cache

    rows = []
    for node_id in eligible:
        node = settings.BEACON_QUERY_NODES[node_id]
        rows.append({
            "node_id": node_id,
            "base_url": node.get("base_url", ""),
            "cache": fresh.get(node_id),
        })

    context = {
        "variant": variant,
        "genome_build": genome_build,
        "rows": rows,
        "skipped": skipped,
    }
    return render(request, "beacon/external_beacons.html", context)
