import re
from snpdb.models import Organization
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample

MIN_3_ALPHA = re.compile(r"[a-zA-Z]{3,}")


@search_receiver(
    search_type=Organization,
    pattern=MIN_3_ALPHA,
    example=SearchExample(
        note="3 or more letter of the organisation's name",
        example="institute"
    )
)
def organization_search(search_input: SearchInputInstance):
    yield Organization.objects.filter(active=True).filter(search_input.q_words())
