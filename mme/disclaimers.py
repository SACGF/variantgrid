""" A connected database's disclaimers and terms must be shown with its query results and
posted on our site. The published baseline (ga4gh/mme-apis/disclaimers) is transcribed into
settings.MME_NODES; anything a node returns live with a response supersedes it.
"""
from django.conf import settings


def node_disclaimer(node_id: str) -> dict:
    """ Static {disclaimer, terms} baseline for a connected node. """
    node = (settings.MME_NODES or {}).get(node_id) or {}
    return {"disclaimer": node.get("disclaimer") or "", "terms": node.get("terms") or ""}


def effective_disclaimer(node_id: str, response_disclaimer: str = "",
                         response_terms: str = "") -> dict:
    """ What to display with a given node's results: live response wins over baseline. """
    static = node_disclaimer(node_id)
    return {
        "disclaimer": response_disclaimer or static["disclaimer"],
        "terms": response_terms or static["terms"],
        "source": "response" if (response_disclaimer or response_terms) else "baseline",
    }


def connected_nodes() -> list[dict]:
    """ Every connected node with its published disclaimer - drives the public page. """
    return [{"node_id": node_id, **node_disclaimer(node_id)}
            for node_id in sorted(settings.MME_NODES or {})]


def mme_response_body(body: dict) -> dict:
    """ Attach OUR disclaimer/terms to an outgoing MME response body, so peers can display
        them with our results the same way we do for theirs. """
    if settings.MME_DISCLAIMER:
        body["disclaimer"] = settings.MME_DISCLAIMER
    if settings.MME_TERMS:
        body["terms"] = settings.MME_TERMS
    return body
