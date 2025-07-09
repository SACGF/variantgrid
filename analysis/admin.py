from django.contrib import admin
from django.contrib.auth.models import User
from guardian.admin import GuardedModelAdmin

from analysis import models
from snpdb.admin_utils import ModelAdminBasics


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


@admin.register(models.Analysis)
class AnalysisAdmin(GuardedModelAdmin):
    list_display = ['id', 'user', 'name', 'description', 'template_type', 'visible', 'created']
    list_filter = (AnalysisUserFilter,)
    search_fields = ('id', 'lab_record_id')


@admin.register(models.AnalysisTemplate)
class AnalysisTemplateAdmin(ModelAdminBasics):
    pass


@admin.register(models.AnalysisTemplateRun)
class AnalysisTemplateRunAdmin(ModelAdminBasics):
    pass


@admin.register(models.AnalysisTemplateVersion)
class AnalysisTemplateVersionAdmin(ModelAdminBasics):
    pass


@admin.register(models.VariantTag)
class VariantTagAdmin(ModelAdminBasics):
    pass


@admin.register(models.VariantTagsImport)
class VariantTagsImportAdmin(ModelAdminBasics):
    pass
