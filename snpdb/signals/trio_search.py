from snpdb.models import Trio
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Trio,
    example=SearchExample(
        note="Provide part of the trio's name"
    )
)
def search_trio(search_input: SearchInputInstance):
    q = search_input.q_words()
    yield Trio.filter_for_user(search_input.user).filter(q, cohort__genome_build__in=search_input.genome_builds)
