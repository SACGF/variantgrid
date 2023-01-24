from django.contrib import admin
from django.contrib.auth.models import User

from snpdb.admin_utils import ModelAdminBasics
from snpdb.models import ProcessingStatus
from . import models
from .models import UploadStep, UploadPipeline, UploadedFile, UploadedVCF


class UploadStepStatusFilter(admin.SimpleListFilter):
    title = 'Status'
    parameter_name = 'status'
    default_value = None

    def lookups(self, request, model_admin):
        return ProcessingStatus.choices

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(status=self.value())
        return queryset


@admin.register(UploadStep)
class UploadStepAdmin(ModelAdminBasics):
    list_display = ('id', 'name', 'upload_pipeline', 'status', 'start_date', 'end_date', 'error_message')
    list_filter = (UploadStepStatusFilter,)
    search_fields = ('id', 'name')

    def mark_timed_out(self, request, queryset):
        user: User = request.user
        us: UploadStep
        for us in queryset:
            us.mark_timed_out(user=user)

    mark_timed_out.short_description = "Mark process as timed out"

    actions = [mark_timed_out]


class UploadStepInline(admin.TabularInline):
    model = UploadStep
    fields = ['id', 'name', 'status']
    show_change_link = True

    def has_add_permission(self, request, obj):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


@admin.register(UploadedFile)
class UploadFileAdmin(ModelAdminBasics):
    list_display = ('name', 'file_type', 'path', 'uploaded_file')


@admin.register(UploadPipeline)
class UploadPipelineAdmin(ModelAdminBasics):
    list_display = ('uploaded_file', 'status')
    inlines = (UploadStepInline, )


@admin.register(UploadedVCF)
class UploadedVCFAdmin(ModelAdminBasics):
    list_display = ('id', 'uploaded_file', 'upload_pipeline', 'vcf')
