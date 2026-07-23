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

from beacon.variant_mapping import variant_to_beacon_query_params, parse_beacon_response
from library.log_utils import AdminNotificationBuilder


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
        headers = {"Accept": "application/json"}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        return headers

    def query_g_variants(self, params: dict, timeout: int) -> dict:
        response = requests.get(
            url=f"{self.base_url}/g_variants",
            params=params,
            headers=self._headers,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()


def query_node(node_id: str, params: dict) -> dict:
    """ Query one node, returning {"node_id", "exists", "count", "error"}.
        A slow/dead node yields an error row, never raises to the caller. """
    result = {"node_id": node_id, "exists": False, "count": None, "error": None}
    try:
        data = BeaconClient(node_id).query_g_variants(params, timeout=settings.BEACON_QUERY_TIMEOUT)
        exists, count = parse_beacon_response(data)
        result["exists"] = exists
        result["count"] = count
    except Exception as e:  # timeout / connection / HTTP error -> "unavailable" for this row
        logging.warning("Beacon outbound: node %s failed: %s", node_id, e)
        result["error"] = str(e)
        nb = AdminNotificationBuilder("Beacon outbound query failed")
        nb.add_markdown(f"Node `{node_id}`: {e}")
        nb.send()
    return result


def query_all_nodes(params: dict) -> list[dict]:
    """ Fan out the same coordinate query across all configured nodes concurrently, each
        bounded by BEACON_QUERY_TIMEOUT, so one slow node never blocks the others (§9.5). """
    node_ids = list(settings.BEACON_QUERY_NODES.keys())
    if not node_ids:
        return []
    with ThreadPoolExecutor(max_workers=len(node_ids)) as executor:
        return list(executor.map(lambda nid: query_node(nid, params), node_ids))


def query_external_beacons_for_variant(variant, genome_build) -> list[dict]:
    """ Build the public coordinate params for a Variant and fan out to all nodes.
        Returns [] for symbolic variants (no explicit coordinate to send). """
    params = variant_to_beacon_query_params(variant, genome_build)
    if params is None:
        return []
    return query_all_nodes(params)
