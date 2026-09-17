from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from pathtests.models import PathologyTest, PathologyTestVersion, get_cases_qs


#@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class PathologyTestAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        active = self.forwarded.get('active', None)
        qs = PathologyTest.objects.all()
        if active:
            qs = qs.filter(activepathologytestversion__isnull=False)
        return qs


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class PathologyTestVersionAutocompleteView(AutocompleteView):
    fields = ['pathology_test__name']

    def get_user_queryset(self, user):
        active = self.forwarded.get('active', None)
        qs = PathologyTestVersion.objects.all()
        if active:
            qs = qs.filter(activepathologytestversion__isnull=False)
        return qs


class CaseAutocompleteView(AutocompleteView):
    """ Matches name or the LIMS code - most cases have no name (see Case.__str__) """
    fields = ['name', 'external_pk__code']

    def get_user_queryset(self, user):
        return get_cases_qs()

    def sort_queryset(self, qs):
        return qs.order_by("-pk")
