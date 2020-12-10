from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis.models import Analysis, AnalysisTemplate, Q, AnalysisTemplateVersion
from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class AnalysisAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = Analysis.filter_for_user(user)
        return qs.filter(template_type__isnull=True)  # Hide templates


class AnalysisTemplateAutocompleteView(AutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        return AnalysisTemplate.filter(user,
                                       requires_sample_somatic=self.forwarded.get("requires_sample_somatic"),
                                       requires_sample_gene_list=self.forwarded.get("requires_sample_gene_list"),
                                       class_name=self.forwarded.get("class_name"),
                                       atv_kwargs={"appears_in_autocomplete": True})
