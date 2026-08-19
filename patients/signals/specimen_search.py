from functools import partial

from patients.models import Extraction, Specimen
from snpdb.search import (
    HAS_3_ANY,
    SearchExample,
    SearchInputInstance,
    SearchResult,
    SearchResultMatchStrength,
    search_receiver,
)


def _named_result(obj, search_string: str) -> SearchResult:
    """ The whole reference typed in is someone naming a row rather than browsing for one, so it can
        be jumped straight to even where a partial match brings back other records """
    match_strength = SearchResultMatchStrength.STRONG_MATCH
    if obj.reference_id and obj.reference_id.casefold() == search_string.casefold():
        match_strength = SearchResultMatchStrength.ID_MATCH
    return SearchResult(obj.preview, match_strength=match_strength)


def _found_via_specimen_result(extraction) -> SearchResult:
    """ Reached through its specimen rather than named itself - still worth showing, but it doesn't
        stand between the specimen reference and the specimen it names """
    return SearchResult(extraction.preview, match_strength=SearchResultMatchStrength.FUZZY_MATCH)


@search_receiver(
    search_type=Specimen,
    # HAS_3_ANY rather than alpha - a TSO 500 specimen reference is all digits
    pattern=HAS_3_ANY,
    example=SearchExample(
        note="Specimen reference",
        examples=["2600000001"]
    )
)
def specimen_search(search_input: SearchInputInstance):
    qs = Specimen.filter_for_user(search_input.user).filter(search_input.q_words('reference_id'))
    yield qs, partial(_named_result, search_string=search_input.search_string.strip())


@search_receiver(
    search_type=Extraction,
    pattern=HAS_3_ANY,
    example=SearchExample(
        note="Extraction reference, or the reference of its specimen",
        examples=["2600000001C"]
    )
)
def extraction_search(search_input: SearchInputInstance):
    """ An extraction may be referred to by its own reference (eg a container suffix), or its specimen's """
    qs = Extraction.filter_for_user(search_input.user)
    named = qs.filter(search_input.q_words('reference_id'))
    yield named, partial(_named_result, search_string=search_input.search_string.strip())

    via_specimen = qs.filter(search_input.q_words('specimen__reference_id')).exclude(pk__in=named)
    yield via_specimen, _found_via_specimen_result
