from functools import partial
from typing import Optional, List
from django.db.models.functions import Lower
from library.preview_request import PreviewProxyModel
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ALPHA_MIN, SearchResult, \
    SearchResultMatchStrength, SearchMessage
import re


def validate_ontology(term: OntologyTerm, preview_proxy: Optional[PreviewProxyModel] = None) -> SearchResult:
    messages: List[SearchMessage] = []
    if term.ontology_service not in OntologyService.LOCAL_ONTOLOGY_PREFIXES:
        messages = [SearchMessage(f"We do not store {term.ontology_service} locally. Ontology page will only provide external links.")]
    elif term.is_stub:
        messages = [SearchMessage("We do not have this term in our database")]
    elif term.status == OntologyTermStatus.DEPRECATED:
        messages = [SearchMessage("This term is obsolete")]
    elif term.status == OntologyTermStatus.NON_CONDITION:
        messages = [SearchMessage("Note this term is not a suitable value for a classification's condition")]
    preview = term.preview
    if preview_proxy:
        preview.category = preview_proxy.preview_category()
        preview.icon = preview_proxy.preview_icon()
    return SearchResult(preview, messages=messages)


ONTOLOGY_TERM_PATTERN = re.compile(r"^(MONDO|OMIM|MIM|HPO|HP|DOID|ORPHANET)\s*:\s*[0-9]+$", re.IGNORECASE)


@search_receiver(
    search_type=OntologyTerm,
    pattern=ONTOLOGY_TERM_PATTERN,
    sub_name="ID",
    example=SearchExample(
        note="Search by the term's identifier, supports MONDO, OMIM, HP",
        examples=["MONDO:0010726", "OMIM:616299", "HP:0001332"]
    ),
    match_strength=SearchResultMatchStrength.ID_MATCH
)
def ontology_search_id(search_input: SearchInputInstance):
    # search by ID
    try:
        term = OntologyTerm.get_or_stub(search_input.search_string)
        if term.ontology_service not in OntologyService.CONDITION_ONTOLOGIES:
            pass
        else:
            yield validate_ontology(term)
    except ValueError:
        # might not be a valid ontology, there's a lot of text that passes the search_input
        pass


@search_receiver(
    search_type=OntologyTerm,
    pattern=re.compile(r"^HGNC\s*:\s*(.*)$"),
    sub_name="Gene Disease Relationships",
    admin_only=True,
    example=SearchExample(
        note="HGNC: followed by a gene symbol",
        examples=["HGNC:BRCA1"]
    ),
    match_strength=SearchResultMatchStrength.STRONG_MATCH
)
def ontology_search_hgnc(search_input: SearchInputInstance):
    try:
        term = OntologyTerm.get_or_stub(search_input.search_string)
        preview = term.preview
        preview.icon = 'fa fa-dna'
        preview.category = "Gene Disease Relationships"
        yield preview
    except ValueError:
        # might not be a valid ontology, there's a lot of text that passes the search_input
        pass


def _ontology_search_name(search_input: SearchInputInstance, ontology_service: OntologyService, preview_proxy: PreviewProxyModel):
    # search by text (but not if matches Gene Symbol - better solution would be to match but give gene symbol higher priority)

    # instead of filtering out if gene symbol exists, return it, but the match strength of name matches is low, and gene symbol matches should be ID
    # if not GeneSymbol.objects.filter(symbol=search_input.search_string).exists():
    qs = OntologyTerm.objects.filter(ontology_service=ontology_service).filter(search_input.q_words()).order_by('status', Lower('name'), 'index')

    # unless we're specifically searching for obsolete, filter them out
    if 'obsolete' not in search_input.search_string:
        qs = qs.exclude(status=OntologyTermStatus.DEPRECATED)

    # limit results to 20 for each kind, need to give the user an overall warning that we're doing this
    yield qs, partial(validate_ontology, preview_proxy=preview_proxy)


# Note that we're quite lucky in that HPO, OMIM, MONDO, Ontology Term are 3 search
OMIM_PREVIEW = PreviewProxyModel("OMIM", OntologyTerm.preview_icon())
MONDO_PREVIEW = PreviewProxyModel("MONDO", OntologyTerm.preview_icon())
HPO_PREVIEW = PreviewProxyModel("HPO", OntologyTerm.preview_icon())


@search_receiver(
    search_type=OMIM_PREVIEW,
    pattern=HAS_3_ALPHA_MIN,
    example=SearchExample(
        note="3 or more letters of the term's name",
        examples=["symphal"]
    ),
    match_strength=SearchResultMatchStrength.FUZZY_MATCH
)
def omim_name_search(search_input: SearchInputInstance):
    return _ontology_search_name(search_input=search_input, ontology_service=OntologyService.OMIM, preview_proxy=OMIM_PREVIEW)


@search_receiver(
    search_type=MONDO_PREVIEW,
    pattern=HAS_3_ALPHA_MIN,
    example=SearchExample(
        note="3 or more letters of the term's name",
        examples=["Rett syndrome"]
    ),
    match_strength=SearchResultMatchStrength.FUZZY_MATCH
)
def mondo_name_search(search_input: SearchInputInstance):
    return _ontology_search_name(search_input=search_input, ontology_service=OntologyService.MONDO, preview_proxy=MONDO_PREVIEW)


@search_receiver(
    search_type=HPO_PREVIEW,
    pattern=HAS_3_ALPHA_MIN,
    example=SearchExample(
        note="3 or more letters of the term's name",
        examples=["pain"]
    ),
    match_strength=SearchResultMatchStrength.FUZZY_MATCH
)
def hpo_name_search(search_input: SearchInputInstance):
    return _ontology_search_name(search_input=search_input, ontology_service=OntologyService.HPO, preview_proxy=HPO_PREVIEW)
