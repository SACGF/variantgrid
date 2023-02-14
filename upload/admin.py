from django.contrib import admin
from django.contrib.auth.models import User

from snpdb.admin_utils import ModelAdminBasics, admin_list_column
from snpdb.models import ProcessingStatus
from snpdb.user_settings_manager import UserSettingsManager
from . import models
from .models import UploadStep, UploadPipeline, UploadedFile, UploadedVCF, UploadedClassificationImport


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
    list_display = ('id', 'uploaded_file_name', 'uploaded_file_date', 'status')
    inlines = (UploadStepInline, )

    @admin_list_column('File Name', order_field='uploaded_file__name')
    def uploaded_file_name(self, obj: UploadPipeline):
        if uploaded_file := obj.uploaded_file:
            return uploaded_file.name

    @admin_list_column('File Date', order_field='uploaded_file__created')
    def uploaded_file_date(self, obj: UploadPipeline):
        if uploaded_file := obj.uploaded_file:
            timezoned = uploaded_file.created.astimezone(UserSettingsManager.get_user_timezone())
            return f"{timezoned.strftime('%Y-%m-%d %H:%M:%S %z')}"


@admin.register(UploadedVCF)
class UploadedVCFAdmin(ModelAdminBasics):
    list_display = ('id', 'uploaded_file', 'upload_pipeline', 'vcf')


@admin.register(UploadedClassificationImport)
class UploadedClassificationImportAdmin(ModelAdminBasics):
    list_display = ('id', 'uploaded_file_name', 'uploaded_file_date', 'classification_import', 'classification_import_allele_info_count')

    @admin_list_column('File Name', order_field='uploaded_file__name')
    def uploaded_file_name(self, obj: UploadedClassificationImport):
        if uploaded_file := obj.uploaded_file:
            return uploaded_file.name

    @admin_list_column('File Date', order_field='uploaded_file__created')
    def uploaded_file_date(self, obj: UploadedClassificationImport):
        if uploaded_file := obj.uploaded_file:
            timezoned = uploaded_file.created.astimezone(UserSettingsManager.get_user_timezone())
            return f"{timezoned.strftime('%Y-%m-%d %H:%M:%S %z')}"

    @admin_list_column('Outstanding Allele Info Count')
    def classification_import_allele_info_count(self, obj: UploadedClassificationImport):
        if classification_import := obj.classification_import:
            return classification_import.importedalleleinfo_set.count()
