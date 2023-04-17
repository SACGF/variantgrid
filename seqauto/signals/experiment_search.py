from django.dispatch import receiver
from seqauto.models import Experiment
from snpdb.search2 import SearchResponse, SearchInput, search_signal


@receiver(search_signal, sender=SearchInput)
def search_experiment(search_string: str, **kwargs) -> SearchResponse:
    search_response = SearchResponse()
    search_response.extend(Experiment.objects.filter(name__icontains=search_string))
    return search_response
