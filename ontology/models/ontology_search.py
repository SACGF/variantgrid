from typing import List, Optional
from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus
from snpdb.search2 import search_receiver, HAS_ALPHA_PATTERN, \
    SearchInputInstance, SearchExample
import re


def validate_ontology(term: OntologyTerm) -> Optional[List[str]]:
    if term.is_stub:
        return ["We do not have this term in our database"]
    elif term.status == OntologyTermStatus.DEPRECATED:
        return ["This term is obsolete"]
    elif term.status == OntologyTermStatus.NON_CONDITION:
        return ["Note this term is not a suitable value for condition"]


ONTOLOGY_TERM_PATTERN = re.compile(r"\w+[:_]\s*.*")


@search_receiver(
    search_type=OntologyTerm,
    pattern=ONTOLOGY_TERM_PATTERN,
    sub_name="Ontology by ID",
    example=SearchExample(
        note="Search by the term's identifier, supports MONDO, OMIM, HP",
        example="MONDO:0010726"
    )
)
def ontology_search_id(search_input: SearchInputInstance):
    # search by ID
    try:
        term = OntologyTerm.get_or_stub(search_input.search_string)
        yield term, validate_ontology(term)
    except ValueError:
        # might not be a valid ontology, there's a lot of text that passes the search_input
        pass


@search_receiver(
    search_type=OntologyTerm,
    pattern=HAS_ALPHA_PATTERN,
    sub_name="Ontology by name",
    example=SearchExample(
        note="Search by part of the term's name",
        example="Rett syndrome"
    )
)
def ontology_search_name(search_input: SearchInputInstance):
    # search by text (but not if matches Gene Symbol - better solution would be to match but give gene symbol higher priority)
    if not GeneSymbol.objects.filter(symbol=search_input.search_string).exists():

        for ontology_service in [OntologyService.MONDO, OntologyService.OMIM, OntologyService.HPO]:

            qs = OntologyTerm.objects.filter(ontology_service=ontology_service).order_by('status', 'name', 'index')
            for word in search_input.search_words:
                qs = qs.filter(name__icontains=word)

            # unless we're specifically searching for obsolete, filter them out
            if 'obsolete' not in search_input.search_string:
                qs = qs.exclude(status=OntologyTermStatus.DEPRECATED)

            # limit results to 20 for each kind, need to give the user an overall warning that we're doing this
            for obj in qs[0:20]:
                yield obj, validate_ontology(obj)
