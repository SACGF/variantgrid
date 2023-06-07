from snpdb.models import Sample
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ANY


@search_receiver(
    search_type=Sample,
    pattern=HAS_3_ANY,
    example=SearchExample(note="3 or more letters/numbers of the sample name")
)
def sample_search(search_input: SearchInputInstance):
    yield Sample.filter_for_user(search_input.user).filter(search_input.q_words())
