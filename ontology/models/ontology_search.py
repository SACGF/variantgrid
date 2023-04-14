from typing import Any, List
from django.dispatch import receiver
from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService
from snpdb.search2 import search_signal, SearchInput, SearchResponse


def validate_ontology(term: OntologyTerm) -> List[str]:
    if term.is_stub:
        return ["We do not have this term in our database"]
    elif term.is_obsolete:
        return ["This term is obsolete"]

@receiver(search_signal, sender=SearchInput)
def search_ontology(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response = SearchResponse()

    valid = False
    # search by ID
    try:
        if search_input.matches_pattern(r"\w+[:_]\s*.*"):
            term = OntologyTerm.get_or_stub(search_input.search_string)
            response.add(term, validate_ontology(term))
    except ValueError:
        # might not be a valid ontology, there's a lot of text that passes the search_input
        pass

    # search by text (but not if matches Gene Symbol - better solution would be to match but give gene symbol higher priority)
    if search_input.matches_has_alpha() and not GeneSymbol.objects.filter(symbol=search_input.search_string).exists():
        valid = True

        qs = OntologyTerm.objects.exclude(ontology_service=OntologyService.HGNC). \
            filter(name__icontains=search_input.search_string). \
            order_by('ontology_service', 'name', 'index')

        # unless we're specifically searching for obsolete, filter them out
        if 'obsolete' not in search_input.search_string:
            qs = qs.exclude(name__icontains='obsolete')

        for obj in qs:
            messages = ["This term is obsolete"] if obj.is_obsolete else None
            response.add(obj, messages=messages)

    return response if valid else None
