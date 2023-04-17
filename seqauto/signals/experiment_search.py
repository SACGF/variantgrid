from typing import Any

from django.dispatch import receiver
from seqauto.models import Experiment
from snpdb.search2 import SearchResponse, SearchInput, search_signal


@receiver(search_signal, sender=SearchInput)
def experiment_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_has_alpha():
        search_response = SearchResponse()
        search_response.extend(Experiment.objects.filter(name__icontains=search_input.search_string))
        return search_response
