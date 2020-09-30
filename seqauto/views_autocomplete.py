from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import WEEK_SECS, MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from seqauto.models import QCColumn, EnrichmentKit


@method_decorator(cache_page(WEEK_SECS), name='dispatch')
class QCColumnAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = QCColumn.objects.all()
        qc_type = self.forwarded.get('qc_type')
        if qc_type:
            qs = qs.filter(qc_type=qc_type)
        return qs


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class EnrichmentKitAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = EnrichmentKit.objects.all()
        # Default is to hide obsolete kits
        show_obsolete = self.forwarded.get('show_obsolete', False)
        if not show_obsolete:
            qs = qs.filter(obsolete=False)
        return qs.order_by("manufacturer", "name", "version")
