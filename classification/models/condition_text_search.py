import requests

from library.constants import MINUTE_SECS
from ontology.models import OntologyService, OntologyTerm


def condition_text_search(search_text: str, row_limit: int = 10) -> list[OntologyTerm]:
    if not search_text or search_text.lower() in {"not set", "not specified", "not provided", "n/a"}:
        # Searching for blank returns everything (29916 records, though you will only get row_limit)
        # This is probably not what you want, so return early without API call
        return []

    response = requests.get(
        'https://api.monarchinitiative.org/v3/api/search', {
            "q": search_text,
            "category": "biolink:Disease",
            "limit": row_limit
        }, timeout=MINUTE_SECS).json()

    results = response.get("items") or []

    terms: list[OntologyTerm] = []
    for result in results:
        term_id = result.get("id") or ""
        if OntologyService.is_unsupported_id(term_id):
            # Monarch can return terms from ontologies we don't support (e.g. MPATH) - skip those.
            # A malformed id still falls through to get_or_stub so a genuine problem isn't hidden.
            continue
        terms.append(OntologyTerm.get_or_stub(term_id))
    return terms
