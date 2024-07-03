import requests

from library.constants import MINUTE_SECS
from ontology.models import OntologyTerm


def condition_text_search(search_text: str, row_limit: int = 10) -> list[OntologyTerm]:
    response = requests.get(
        'https://api.monarchinitiative.org/v3/api/search', {
            "q": search_text,
            "category": "biolink:Disease",
            "limit": row_limit
        }, timeout=MINUTE_SECS).json()

    results = response.get("items")

    return [OntologyTerm.get_or_stub(result.get("id")) for result in results]
