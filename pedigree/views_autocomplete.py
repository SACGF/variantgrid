from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from library.constants import MINUTE_SECS
from pedigree.models import Pedigree
from snpdb.views.views_autocomplete import GenomeBuildAutocompleteView


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class PedigreeAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = Pedigree.filter_for_user(user, success_status_only=True)
        return self.filter_to_genome_build(qs, "cohort__genome_build")
