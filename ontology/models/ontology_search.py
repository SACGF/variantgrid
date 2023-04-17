from typing import Any, List
from django.dispatch import receiver
from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus
from snpdb.search2 import search_signal, SearchInput, SearchResponse


def validate_ontology(term: OntologyTerm) -> List[str]:
    if term.is_stub:
        return ["We do not have this term in our database"]
    elif term.status == OntologyTermStatus.DEPRECATED:
        return ["This term is obsolete"]
    elif term.status == OntologyTermStatus.NON_CONDITION:
        return ["Note this term is not a suitable value for condition"]


@receiver(search_signal, sender=SearchInput)
def search_ontology(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response = SearchResponse()
    # search by ID
    try:
        if search_input.matches_pattern(r"\w+[:_]\s*.*"):
            response.add_search_category(OntologyTerm)
            term = OntologyTerm.get_or_stub(search_input.search_string)
            response.add(term, validate_ontology(term))
    except ValueError:
        # might not be a valid ontology, there's a lot of text that passes the search_input
        pass

    # search by text (but not if matches Gene Symbol - better solution would be to match but give gene symbol higher priority)
    if search_input.matches_has_alpha() and not GeneSymbol.objects.filter(symbol=search_input.search_string).exists():
        response.add_search_category(OntologyTerm)

        for ontology_service in [OntologyService.MONDO, OntologyService.OMIM, OntologyService.HPO]:

            qs = OntologyTerm.objects.filter(ontology_service=ontology_service).order_by('status', 'name', 'index')
            for word in search_input.search_words:
                qs = qs.filter(name__icontains=word)

            # unless we're specifically searching for obsolete, filter them out
            if 'obsolete' not in search_input.search_string:
                qs = qs.exclude(status=OntologyTermStatus.DEPRECATED)

            # limit results to 20 for each kind, need to give the user an overall warning that we're doing this
            for obj in qs[0:20]:
                response.add(obj, messages=validate_ontology(obj))

    return response
