from typing import Any, Optional, List

from django.dispatch import receiver

from annotation.models.models_citations import CitationSource, CitationIdNormalized, Citation
from snpdb.search2 import SearchResponseRecordAbstract, SearchInput, search_signal, SearchResponse


class SearchResponseCitation(SearchResponseRecordAbstract[Citation]):

    @classmethod
    def search_type(cls) -> str:
        return "Citation"

    @property
    def messages(self) -> Optional[List[str]]:
        input_string = self.search_input.search_string
        tidy_input = input_string.replace(' ', '').upper()
        if ':' in tidy_input:
            colon_index = tidy_input.index(':')
            search_prefix = CitationSource.from_legacy_code(tidy_input[:colon_index])
            suffix = tidy_input[colon_index+1:]
            if self.record.source != search_prefix or self.record.index != suffix:
                return [f'Normalising "{input_string}" to "{self.record.id}"']
        return None


@receiver(search_signal, sender=SearchInput)
def search_citations(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response: SearchResponse[Citation] = SearchResponse(SearchResponseCitation, search_input)

    try:
        normal_id = CitationIdNormalized.normalize_id(search_input.search_string)
        response.mark_valid_search()

        citation = normal_id.get_or_create()
        response.add(citation)
    except ValueError:
        pass

    return response
