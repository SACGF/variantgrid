""" Build the standard Beacon v2 response envelope and clamp granularity by auth.

Centralising envelope construction + the granularity clamp here means no view can leak
more detail than the requester's tier allows (§7).
"""
from django.conf import settings

# Beacon granularity tiers, lowest disclosure first.
BOOLEAN = "boolean"
COUNT = "count"
RECORD = "record"
GRANULARITY_ORDER = {BOOLEAN: 0, COUNT: 1, RECORD: 2}
GRANULARITIES = (BOOLEAN, COUNT, RECORD)


def _rank(granularity: str) -> int:
    return GRANULARITY_ORDER.get(granularity, 0)


def _min_granularity(a: str, b: str) -> str:
    return a if _rank(a) <= _rank(b) else b


def clamp_granularity(requested: str, authenticated: bool) -> str:
    """ Resolve the granularity to return: min(requested-or-default, allowed-by-auth, max).

        Data scope is enforced separately (filter_for_user -> anonymous sees only the
        public group), so the granularity ceiling is currently the configured max for both
        anonymous and authenticated requesters (open question #4 - no OIDC scoping yet).
        `authenticated` is threaded through so a future security pass can lower the anon
        ceiling here without touching call sites. """
    config = settings.BEACON_CONFIG
    ceiling = config.get("max_granularity", RECORD)
    chosen = requested or config.get("default_granularity", BOOLEAN)
    if chosen not in GRANULARITY_ORDER:
        chosen = config.get("default_granularity", BOOLEAN)
    allowed_by_auth = ceiling  # both tiers currently allowed up to max; see docstring
    return _min_granularity(_min_granularity(chosen, allowed_by_auth), ceiling)


def beacon_meta(returned_granularity: str, received_request_summary: dict) -> dict:
    config = settings.BEACON_CONFIG
    return {
        "beaconId": config["beacon_id"],
        "apiVersion": config["api_version"],
        "returnedGranularity": returned_granularity,
        "receivedRequestSummary": received_request_summary,
    }


def query_response(returned_granularity: str, received_request_summary: dict, *,
                   exists: bool, num_total_results: int, result_sets: list[dict]) -> dict:
    """ Assemble a g_variants query response, trimming detail to the returned granularity:
        - boolean: meta + responseSummary.exists only
        - count:   + responseSummary.numTotalResults
        - record:  + response.resultSets[] (one per dataset)
    """
    response = {
        "meta": beacon_meta(returned_granularity, received_request_summary),
        "responseSummary": {"exists": exists},
    }
    if _rank(returned_granularity) >= _rank(COUNT):
        response["responseSummary"]["numTotalResults"] = num_total_results
    if _rank(returned_granularity) >= _rank(RECORD):
        response["response"] = {"resultSets": result_sets}
    return response


def info_response(payload: dict) -> dict:
    """ Framework (non-query) endpoints wrap their payload in a minimal meta. """
    config = settings.BEACON_CONFIG
    return {
        "meta": {
            "beaconId": config["beacon_id"],
            "apiVersion": config["api_version"],
            "returnedSchemas": [],
        },
        "response": payload,
    }
