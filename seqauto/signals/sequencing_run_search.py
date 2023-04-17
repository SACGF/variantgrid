from typing import Any, Iterable

from django.dispatch import receiver

from seqauto.models import SequencingRun
from snpdb.search2 import SearchInput, SearchResponse, search_signal

SEQUENCING_RUN_REGEX = r"\d{6}[_-](NS|NB|M|D|SN|K|ST|A)(.{3,7})_\d{4}_(0{9}-.{5}|.{10})"


@receiver(search_signal, sender=SearchInput)
def sequencing_run_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_pattern(SEQUENCING_RUN_REGEX):
        # TODO test
        response = SearchResponse(SequencingRun)
        response.extend(SequencingRun.objects.filter(name__icontains=search_input.search_string))
        return response
