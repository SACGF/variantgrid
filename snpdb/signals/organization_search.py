import re
from typing import Type, Any

from django.dispatch import receiver

from snpdb.models import Organization
from snpdb.search2 import search_signal, SearchInput, SearchResponse

MIN_3_ALPHA = re.compile(r"[a-zA-Z]{3,}")



@receiver(search_signal, sender=SearchInput)
def organization_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_pattern(MIN_3_ALPHA):
        response = SearchResponse(Organization)
        response.extend(Organization.objects.filter(active=True).filter(name__icontains=search_input.search_string.upper()))
        return response
