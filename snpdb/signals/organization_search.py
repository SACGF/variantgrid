from snpdb.models import Organization
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ALPHA_MIN


@search_receiver(
    search_type=Organization,
    pattern=HAS_3_ALPHA_MIN,
    example=SearchExample(
        note="3 or more letters of the organisation's name",
        examples=["institute"]
    )
)
def organization_search(search_input: SearchInputInstance):
    yield Organization.objects.filter(active=True).filter(search_input.q_words())
