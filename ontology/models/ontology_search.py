from typing import Any, Optional
from django.dispatch import receiver
from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService
from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse


class SearchResponseOntology(SearchResponseRecordAbstract[OntologyTerm]):

    @classmethod
    def search_type(cls) -> str:
        return "Ontology"


@receiver(search_signal, sender=SearchInput)
def search_ontology(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response: SearchResponse[OntologyTerm] = SearchResponse(SearchResponseOntology)

    if search_input.matches_pattern(r"\w+:\s*.*"):
        response.mark_valid_search()

        term = OntologyTerm.get_or_stub(search_input.search_string)
        if not term.is_stub:
            response.add(term)

    if search_input.matches_has_alpha() and not GeneSymbol.objects.filter(symbol=search_input.search_string).exists():
        qs = OntologyTerm.objects.exclude(ontology_service=OntologyService.HGNC). \
            filter(name__icontains=search_input.search_string). \
            order_by('ontology_service', 'name', 'index')

        # unless we're specifically searching for obsolete, filter them out
        if 'obsolete' not in search_input.search_string:
            qs = qs.exclude(name__icontains='obsolete')

        response.extend(qs)

    return response
