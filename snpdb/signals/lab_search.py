
from snpdb.models import Lab
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample
import re


MIN_3_ALPHA = re.compile(r"[a-zA-Z]{3,}")


@search_receiver(
    search_type=Lab,
    pattern=MIN_3_ALPHA,
    example=SearchExample(
        note="3 or more letter of the lab's name",
        examples=["molecular"]
    )
)
def lab_search(search_input: SearchInputInstance):
    yield Lab.objects.filter(organization__active=True).filter(search_input.q_words())

