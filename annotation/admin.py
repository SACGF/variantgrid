from datetime import timedelta

from django.contrib import admin
from django.db.models import QuerySet

from annotation import models
from annotation.models import Citation, CitationFetchRequest
from snpdb.admin_utils import ModelAdminBasics, admin_action

admin.site.register(models.AnnotationRun)
admin.site.register(models.AnnotationVersion)
admin.site.register(models.ClinVar)
admin.site.register(models.VariantAnnotationVersion)

class HasErrorFilter(admin.SimpleListFilter):
    title = "Has Errors"
    parameter_name = "error_mode"

    def lookups(self, request, model_admin):
        return [("has_error", "Has Errors")]

    def queryset(self, request, queryset: QuerySet[Citation]):
        if self.value() == "has_error":
            queryset = queryset.filter(error__isnull=False)
        return queryset


@admin.register(Citation)
class CitationAdmin(ModelAdminBasics):
    list_display = ('id', 'last_loaded', 'error', 'old_id')
    list_filter = ('source', HasErrorFilter)
    search_fields = ('id', 'source')

    @admin_action("Force Refresh")
    def force_refresh(self, request, queryset: QuerySet[Citation]):
        CitationFetchRequest.fetch_all_now(queryset, cache_age=timedelta(seconds=0))
