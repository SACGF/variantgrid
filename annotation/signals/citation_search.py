from typing import Any
from django.dispatch import receiver
from annotation.models import CitationFetchRequest
from annotation.models.models_citations import CitationSource, CitationIdNormalized
from snpdb.search2 import SearchInput, search_signal, SearchResponse


@receiver(search_signal, sender=SearchInput)
def search_citations(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    try:
        normal_id = CitationIdNormalized.normalize_id(search_input.search_string)
        response = SearchResponse("Citation")

        citation = normal_id.get_or_create()
        CitationFetchRequest.fetch_all_now([normal_id])

        messages = []
        input_string = search_input.search_string
        tidy_input = input_string.replace(' ', '').upper()
        if ':' in tidy_input:
            colon_index = tidy_input.index(':')
            search_prefix = CitationSource.from_legacy_code(tidy_input[:colon_index])
            suffix = tidy_input[colon_index + 1:]
            if citation.source != search_prefix or citation.index != suffix:
                messages.append(f'Normalising "{input_string}" to "{citation.id}"')
        if citation.error:
            messages.append("Could not retrieve citation")

        response.add(citation, messages=messages)
        return response
    except ValueError:
        pass

