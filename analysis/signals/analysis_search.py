from typing import Any
from django.dispatch import receiver
from analysis.models import Analysis
from snpdb.search2 import search_signal, SearchInput, SearchResponse
import re


ANALYSIS_PREFIX_PATTERN = re.compile(r"^a(\d+)$")


@receiver(search_signal, sender=SearchInput)
def search_analysis(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if m := search_input.matches_pattern(ANALYSIS_PREFIX_PATTERN):
        response = SearchResponse(Analysis)
        response.extend(Analysis.objects.filter(pk=m.group(1)))

        return response
