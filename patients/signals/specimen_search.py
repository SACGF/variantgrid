from patients.models import Extraction, Specimen
from snpdb.search import HAS_3_ANY, SearchExample, SearchInputInstance, search_receiver


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
    yield Specimen.filter_for_user(search_input.user).filter(search_input.q_words('reference_id'))


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
    q_reference = search_input.q_words('reference_id') | search_input.q_words('specimen__reference_id')
    yield Extraction.filter_for_user(search_input.user).filter(q_reference)
