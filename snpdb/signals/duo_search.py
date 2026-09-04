from snpdb.models import Duo
from snpdb.search import SearchExample, SearchInputInstance, search_receiver


@search_receiver(
    search_type=Duo,
    example=SearchExample(
        note="Provide part of the duo's name"
    )
)
def search_duo(search_input: SearchInputInstance):
    q = search_input.q_words()
    yield Duo.filter_for_user(search_input.user).filter(q, cohort__genome_build__in=search_input.genome_builds)
