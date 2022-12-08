from datetime import timedelta

from django.contrib import admin
from django.db.models import QuerySet

from annotation import models
from annotation.models import Citation2, CitationFetchRequest
from snpdb.admin_utils import ModelAdminBasics, admin_action

admin.site.register(models.AnnotationRun)
admin.site.register(models.AnnotationVersion)
admin.site.register(models.ClinVar)
admin.site.register(models.VariantAnnotationVersion)


@admin.register(Citation2)
class Citation2Admin(ModelAdminBasics):
    list_display = ('id', 'last_loaded', 'error', 'old_id')
    list_filter = ('source',)
    search_fields = ('id', 'source')

    @admin_action("Force Refresh")
    def force_refresh(self, request, queryset: QuerySet[Citation2]):
        CitationFetchRequest.fetch_all_now(queryset, cache_age=timedelta(seconds=0))
