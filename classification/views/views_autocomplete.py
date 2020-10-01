from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from classification.models import EvidenceKey


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class EvidenceKeyAutocompleteView(AutocompleteView):
    fields = ['key']

    def get_user_queryset(self, user):
        return EvidenceKey.objects.all()
