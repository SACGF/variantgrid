from snpdb.models import Quad
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Quad,
    example=SearchExample(
        note="Provide part of the quad's name"
    )
)
def search_quad(search_input: SearchInputInstance):
    q = search_input.q_words()
    yield Quad.filter_for_user(search_input.user).filter(q, cohort__genome_build__in=search_input.genome_builds)
