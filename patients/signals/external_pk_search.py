from django.core.exceptions import ObjectDoesNotExist
from rest_framework.exceptions import PermissionDenied

from patients.models import ExternalPK
from snpdb.search2 import search_receiver, HAS_ALPHA_PATTERN, \
    SearchInputInstance, SearchExample


@search_receiver(
    search_type=ExternalPK, # FIXME, not really appropriate to call this ExternalPK
    pattern=HAS_ALPHA_PATTERN,
    example=SearchExample(
        note="Search on HelixID or SAPOrderNumber"
    )
)
def search_external_pk(search_input: SearchInputInstance):
    # TODO test me
    # Returns related objects
    RELATED_OBJECT_FIELDS = ["case", "pathologytestorder", "patient"]
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
                yield obj
            except ObjectDoesNotExist:
                pass