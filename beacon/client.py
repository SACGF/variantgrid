""" Outbound Beacon client (§9).

Query external Beacons for a variant we're viewing: we send only a public coordinate
(referenceName/start/referenceBases/alternateBases/assemblyId) and read back
presence/count. A live read-enrichment, not a submission of our data - so unlike MME
outbound there is no draft/confirm workflow and no persisted candidate rows, just a
cached (variant, node) presence/count (§9.5).

Transport mirrors mme/client.py: a small `requests` wrapper, one instance per node.
"""
import logging
from concurrent.futures import ThreadPoolExecutor

import requests
from django.conf import settings

from beacon.query_targets import eligible_queries
from beacon.variant_mapping import parse_beacon_response

# Identify ourselves to remote Beacons - some hosting WAFs block default client user-agents
# (e.g. the progenetix family 403s Python clients). The Beacon spec itself requires no UA.
USER_AGENT = "VariantGrid-Beacon/1.0 (+https://variantgrid.com)"


class BeaconClient:
    """ Minimal external-Beacon client. One instance per configured node. """

    def __init__(self, node_id: str):
        node = settings.BEACON_QUERY_NODES[node_id]
        self.node_id = node_id
        self.base_url = node["base_url"].rstrip("/")
        self.api_version = node.get("api_version", "v2.0.0")
        self.token = node.get("token")  # optional: open tiers need none

    @property
    def _headers(self) -> dict:
        headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        return headers

    def query_g_variants(self, params: dict, timeout: int) -> dict:
        # We only ever render presence/count, so ask for count granularity - a beacon clamps
        # it down to boolean if that's all its open tier allows, and it avoids pulling records.
        response = requests.get(
            url=f"{self.base_url}/g_variants",
            params={**params, "requestedGranularity": "count"},
            headers=self._headers,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()


def query_node(node_id: str, params: dict) -> dict:
    """ Query one node, returning {"node_id", "exists", "count", "error"}.
        A slow/dead node yields an error row, never raises to the caller. `exists` stays None
        unless the node actually answered presence/absence. """
    result = {"node_id": node_id, "exists": None, "count": None, "error": None}
    try:
        data = BeaconClient(node_id).query_g_variants(params, timeout=settings.BEACON_QUERY_TIMEOUT)
        exists, count = parse_beacon_response(data)
        result["exists"] = exists
        result["count"] = count
    except Exception as e:  # timeout / connection / HTTP error -> "unavailable" for this row
        # A slow/dead external Beacon is an expected condition (§9.5): log it and surface an
        # error row - it is not an admin-worthy incident, so no notification.
        logging.warning("Beacon outbound: node %s failed: %s", node_id, e)
        result["error"] = str(e)
    return result


def query_nodes(node_params: dict[str, dict]) -> list[dict]:
    """ Fan out per-node queries concurrently - each node with its own params, each bounded by
        BEACON_QUERY_TIMEOUT, so one slow node never blocks the others (§9.5). """
    if not node_params:
        return []
    items = list(node_params.items())
    with ThreadPoolExecutor(max_workers=len(items)) as executor:
        return list(executor.map(lambda item: query_node(item[0], item[1]), items))


def query_external_beacons_for_variant(variant, genome_build) -> list[dict]:
    """ Gate the variant against each configured node's target and fan out only to the servers
        whose domain it matches (§9): symbolic CNVs -> copy-number Beacons, SNVs -> sequence
        Beacons. A variant that matches no configured server queries nothing. """
    node_params = eligible_queries(variant, genome_build, settings.BEACON_QUERY_NODES)
    return query_nodes(node_params)
