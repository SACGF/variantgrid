from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_SEARCH_REGEX
from seqauto.models import SequencingRun
from snpdb.search import SearchExample, SearchInputInstance, search_receiver


@search_receiver(
    search_type=SequencingRun,
    pattern=SEQUENCING_RUN_SEARCH_REGEX,
    example=SearchExample(
        "Part of the name of a sequencing run"
    )
)
def sequencing_run_search(search_input: SearchInputInstance):
    yield SequencingRun.objects.filter(search_input.q_words())
