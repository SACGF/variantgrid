from django.contrib import admin
from django.contrib.auth.models import User
from guardian.admin import GuardedModelAdmin

from analysis import models


class AnalysisUserFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'User Filter'
    parameter_name = 'user'
    default_value = None

    def lookups(self, request, model_admin):
        return [(user.id, user.username) for user in User.objects.all()]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(user=self.value())
        return queryset


class AnalysisAdmin(GuardedModelAdmin):
    list_display = ['id', 'user', 'name', 'description', 'template_type', 'visible', 'created']
    list_filter = (AnalysisUserFilter,)
    search_fields = ('id', 'lab_record_id')


admin.site.register(models.Analysis, AnalysisAdmin)
admin.site.register(models.AnalysisTemplate)
admin.site.register(models.AnalysisTemplateRun)
admin.site.register(models.AnalysisTemplateVersion)
admin.site.register(models.VariantTag)
admin.site.register(models.VariantTagsImport)
