import requests
from requests import RequestException

from library.constants import MINUTE_SECS
from ontology.models import OntologyTerm


def condition_text_search(search_text: str, row_limit: int = 10) -> list[OntologyTerm]:
    if not search_text or search_text.lower() in {"not set", "not specified", "not provided", "n/a"}:
        # Searching for blank returns everything (29916 records, though you will only get row_limit)
        # This is probably not what you want, so return early without API call
        return []

    try:
        http_response = requests.get(
            'https://api.monarchinitiative.org/v3/api/search', {
                "q": search_text,
                "category": "biolink:Disease",
                "limit": row_limit
            }, timeout=MINUTE_SECS)
        http_response.raise_for_status()
        response = http_response.json()
    except (RequestException, ValueError):
        # The Monarch API is external and intermittently returns errors or non-JSON error pages
        # (e.g. a 502 HTML page), which makes .json() raise. Treat any such failure as "no results"
        # rather than aborting the whole condition search.
        return []

    results = response.get("items") or []

    terms: list[OntologyTerm] = []
    for result in results:
        try:
            terms.append(OntologyTerm.get_or_stub(result.get("id")))
        except ValueError:
            # The Monarch search can return terms from ontologies we don't support
            # (e.g. MPATH), whose prefix isn't a member of OntologyService - skip those
            continue
    return terms
