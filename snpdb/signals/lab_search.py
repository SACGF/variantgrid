
from snpdb.models import Lab
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ALPHA_MIN


@search_receiver(
    search_type=Lab,
    pattern=HAS_3_ALPHA_MIN,
    example=SearchExample(
        note="3 or more letters of the lab's name",
        examples=["molecular"]
    )
)
def lab_search(search_input: SearchInputInstance):
    yield Lab.objects.filter(organization__active=True).filter(search_input.q_words())
