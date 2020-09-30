from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from annotation.models.models_mim_hpo import MIMMorbidAlias, HPOSynonym
from library.constants import WEEK_SECS
from library.django_utils.autocomplete_utils import AutocompleteView


@method_decorator(cache_page(WEEK_SECS), name='dispatch')
class MIMMorbidAliasAutocompleteView(AutocompleteView):
    fields = ['description']

    def get_user_queryset(self, user):
        return MIMMorbidAlias.objects.all()


@method_decorator(cache_page(WEEK_SECS), name='dispatch')
class HPOSynonymAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        return HPOSynonym.objects.all()
