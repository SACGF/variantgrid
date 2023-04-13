from typing import Any, Optional, List, Union, Type

from django.dispatch import receiver
from django.utils.safestring import SafeString

from annotation.models.models_citations import CitationSource, CitationIdNormalized, Citation
from ontology.models import OntologyTerm
from snpdb.search2 import SearchResponseRecordAbstract, SearchInput, search_signal, SearchResponse


class SearchResponseCitation(SearchResponseRecordAbstract[Citation]):

    @classmethod
    def result_class(cls) -> Type:
        return OntologyTerm

    @property
    def messages(self) -> Optional[List[str]]:
        issues = []
        input_string = self.search_input.search_string
        tidy_input = input_string.replace(' ', '').upper()
        if ':' in tidy_input:
            colon_index = tidy_input.index(':')
            search_prefix = CitationSource.from_legacy_code(tidy_input[:colon_index])
            suffix = tidy_input[colon_index+1:]
            if self.record.source != search_prefix or self.record.index != suffix:
                issues.append(f'Normalising "{input_string}" to "{self.record.id}"')
        if self.record.error:
            issues.append("Could not retrieve summary for citation")

        if issues:
            return issues

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
