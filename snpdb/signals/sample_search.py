from snpdb.models import Sample
from snpdb.search import search_receiver, SearchInputInstance, HAS_ALPHA_PATTERN, SearchExample


@search_receiver(
    search_type=Sample,
    pattern=HAS_ALPHA_PATTERN,
    example=SearchExample(note="Provide part of the sample's name")
)
def sample_search(search_input: SearchInputInstance):
    yield Sample.filter_for_user(search_input.user).filter(search_input.q_words())