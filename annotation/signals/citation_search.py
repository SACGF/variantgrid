from typing import Any

from django.dispatch import receiver

from annotation.models.models_citations import CitationSource, CitationIdNormalized, Citation
from snpdb.search2 import SearchResponseRecordAbstract, SearchInput, search_signal, SearchResponse


class SearchResponseCitation(SearchResponseRecordAbstract[CitationSource]):

    @classmethod
    def search_type(cls) -> str:
        return "Citation"


@receiver(search_signal, sender=SearchInput)
def citation_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response: SearchResponse[Citation] = SearchResponse(SearchResponseCitation)

    try:
        normal_id = CitationIdNormalized.normalize_id(search_input.search_string)
        response.mark_valid_search()

        citation = normal_id.get_or_create()
        response.add(citation)
    except ValueError:
        pass

    return response
