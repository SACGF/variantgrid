from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis.models import Analysis, AnalysisTemplateVersion
from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class AnalysisAutocompleteView(AutocompleteView):
    fields = ['pk', 'name']

    def get_user_queryset(self, user):
        qs = Analysis.filter_for_user(user)
        if template_type := self.forwarded.get('template_type', None):
            qs = qs.filter(template_type=template_type)
        else:
            qs = qs.filter(template_type__isnull=True)  # Hide templates
        return qs


class AnalysisTemplateAutocompleteView(AutocompleteView):
    fields = ['template__name']

    def sort_queryset(self, qs):
        return qs.order_by("template__name", "-version")

    def get_user_queryset(self, user):
        return AnalysisTemplateVersion.filter_for_user(
            user,
            requires_sample_somatic=self.forwarded.get("requires_sample_somatic"),
            requires_sample_gene_list=self.forwarded.get("requires_sample_gene_list"),
            class_name=self.forwarded.get("class_name"),
            appears_in_autocomplete=True)
