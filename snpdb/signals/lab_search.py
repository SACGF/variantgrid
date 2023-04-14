from typing import Type, Any
from django.dispatch import receiver
from snpdb.models import Lab
from snpdb.search2 import search_signal, SearchInput, SearchResponse
import re


MIN_3_ALPHA = re.compile(r"[a-zA-Z]{3,}")



@receiver(search_signal, sender=SearchInput)
def lab_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_pattern(MIN_3_ALPHA):
        response = SearchResponse("Lab")
        response.extend(Lab.objects.filter(organization__active=True).filter(name__icontains=search_input.search_string.upper()))
        return response
