from typing import Any
from django.core.exceptions import ObjectDoesNotExist
from django.dispatch import receiver
from rest_framework.exceptions import PermissionDenied
from patients.models import ExternalPK
from snpdb.search2 import search_signal, SearchInput, SearchResponse


@receiver(search_signal, sender=SearchInput)
def search_external_pk(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    # TODO test me
    # Returns related objects
    if search_input.matches_has_alpha():
        RELATED_OBJECT_FIELDS = ["case", "pathologytestorder", "patient"]
        response = SearchResponse(ExternalPK)
        for external_pk in ExternalPK.objects.filter(code__iexact=search_input.search_string):
            for f in RELATED_OBJECT_FIELDS:
                try:
                    obj = getattr(external_pk, f)
                    try:
                        obj.check_can_view(search_input.user)
                    except PermissionDenied:
                        continue  # Don't add to results
                    except AttributeError:
                        pass  # No permissions - ok to add to results
                    response.append(obj)
                except ObjectDoesNotExist:
                    pass
        return response