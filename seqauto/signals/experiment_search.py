from seqauto.models import Experiment
from snpdb.search import search_receiver, HAS_ALPHA_PATTERN, SearchInputInstance


@search_receiver(search_type=Experiment, pattern=HAS_ALPHA_PATTERN)
def experiment_search(search_input: SearchInputInstance):
    yield Experiment.objects.filter(search_input.q_words())
