from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from library.constants import MINUTE_SECS
from ontology.models import OntologyService, OntologyTerm

MONARCH_SEARCH_URL = 'https://api.monarchinitiative.org/v3/api/search'

# Monarch intermittently returns a gateway error, so retry once before letting the failure through
# (raise_on_status=False so the caller's raise_for_status reports the actual status, not a RetryError)
_retry = Retry(total=1, status_forcelist=[500, 502, 503, 504], allowed_methods=["GET"],
               backoff_factor=1, raise_on_status=False)
_session = Session()
_session.mount("https://", HTTPAdapter(max_retries=_retry))


def condition_text_search(search_text: str, row_limit: int = 10) -> list[OntologyTerm]:
    if not search_text or search_text.lower() in {"not set", "not specified", "not provided", "n/a"}:
        # Searching for blank returns everything (29916 records, though you will only get row_limit)
        # This is probably not what you want, so return early without API call
        return []

    http_response = _session.get(
        MONARCH_SEARCH_URL, {
            "q": search_text,
            "category": "biolink:Disease",
            "limit": row_limit
        }, timeout=MINUTE_SECS)
    http_response.raise_for_status()
    try:
        response = http_response.json()
    except ValueError as ve:
        # Monarch can return 200 with a non-JSON body (e.g. a proxy error page) - say who sent what,
        # otherwise this propagates as a bare JSONDecodeError that looks like a bug in our parsing
        content_type = http_response.headers.get("Content-Type")
        raise ValueError(f"Monarch search returned non-JSON ({content_type}): "
                         f"{http_response.text[:200]}") from ve

    results = response.get("items") or []

    terms: list[OntologyTerm] = []
    for result in results:
        term_id = result.get("id") or ""
        if not OntologyService.is_supported_id(term_id):
            # Monarch can return ids from ontologies we don't support (e.g. MPATH) or the odd
            # malformed id - skip rather than letting one bad result abort the whole search.
            continue
        terms.append(OntologyTerm.get_or_stub(term_id))
    return terms
