import re

from seqauto.models import SequencingRun
from snpdb.search import search_receiver, SearchInputInstance, SearchExample

SEQUENCING_RUN_REGEX = re.compile(r"\d{6}[_-](NS|NB|M|D|SN|K|ST|A)(.{3,7})_\d{4}_(0{9}-.{5}|.{10})")


@search_receiver(
    search_type=SequencingRun,
    pattern=SEQUENCING_RUN_REGEX,
    example=SearchExample(
        "Part of the name of a sequencing run"
    )
)
def sequencing_run_search(search_input: SearchInputInstance):
    yield SequencingRun.objects.filter(search_input.q_words())
