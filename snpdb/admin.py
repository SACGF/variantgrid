import csv

from django.contrib import admin
from django.db.models import AutoField, ForeignKey
from django.http import HttpResponse
from guardian.admin import GuardedModelAdmin

from snpdb import models
from snpdb.models import Allele
from snpdb.models.models_genome import GenomeBuild
from snpdb.models_admin_forms import LabAdmin, OrganizationAdmin


class AdminExportCsvMixin:
    def export_as_csv(self, request, queryset):

        meta = self.model._meta
        field_names = [field.name for field in meta.fields]

        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={}.csv'.format(meta)
        writer = csv.writer(response)

        writer.writerow(field_names)
        for obj in queryset:
            row = writer.writerow([getattr(obj, field) for field in field_names])
        return response

    export_as_csv.short_description = "Export Selected"

    def _get_readonly_fields(self, request, obj=None):
        readonly = []
        for f in self.model._meta.fields:
            if isinstance(f, (AutoField, ForeignKey)):
                readonly.append(f.name)
        return readonly

    def _get_fields(self, request, obj=None, **kwargs):
        readonly = []
        mutable = []
        for f in self.model._meta.fields:
            if isinstance(f, (AutoField, ForeignKey)):
                readonly.append(f.name)
            else:
                mutable.append(f.name)

        return readonly + mutable


class ModelAdminBasics(admin.ModelAdmin, AdminExportCsvMixin):
    # wanted to call this BaseModelAdmin but that was already taken
    actions = ["export_as_csv"]

    def get_readonly_fields(self, request, obj=None):
        return self._get_readonly_fields(request=request, obj=obj)

    def get_fields(self, request, obj=None):
        return self._get_fields(request=request, obj=obj)


class GuardedModelAdminBasics(GuardedModelAdmin, AdminExportCsvMixin):
    actions = ["export_as_csv"]

    def get_readonly_fields(self, request, obj=None):
        return self._get_readonly_fields(request=request, obj=obj)

    def get_fields(self, request, obj=None):
        return self._get_fields(request=request, obj=obj)


class AlleleAdmin(admin.ModelAdmin, AdminExportCsvMixin):

    def validate(self, request, queryset):
        allele: Allele
        for allele in queryset:
            allele.validate()

    validate.short_description = 'Validate'
    actions = ["export_as_csv", validate]
    search_fields = ('id', 'clingen_allele__id')
    # flag collection and clingen allele blow up as drop downs as there's so many choices in them
    fieldsets = (
        ('No Fields', {'fields': ()}),
    )


class DefaultBuildFilter(admin.SimpleListFilter):
    title = 'Default Genome Build Filter'
    parameter_name = 'default_genome_build'
    default_value = None

    def lookups(self, request, model_admin):
        return [(gb.name, gb.name) for gb in GenomeBuild.objects.filter(enabled=True).order_by('name')]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(default_genome_build=self.value())
        return queryset


class UserSettingsOverrideAdmin(admin.ModelAdmin, AdminExportCsvMixin):
    list_display = ('user', 'default_lab', 'default_genome_build', 'email_weekly_updates', 'email_discordance_updates')
    list_filter = (DefaultBuildFilter,)
    search_fields = ('user__username', 'user__first_name', 'user__last_name')
    actions = ["export_as_csv"]


admin.site.register(models.Allele, AlleleAdmin)
admin.site.register(models.CachedGeneratedFile, ModelAdminBasics)
admin.site.register(models.Cohort, ModelAdminBasics)
admin.site.register(models.CohortGenotypeCollection, ModelAdminBasics)
admin.site.register(models.CohortSample, ModelAdminBasics)
admin.site.register(models.CustomColumn, ModelAdminBasics)
admin.site.register(models.CustomColumnsCollection, ModelAdminBasics)
admin.site.register(models.GenomeBuild, ModelAdminBasics)
admin.site.register(models.GenomicIntervalsCollection, ModelAdminBasics)
admin.site.register(models.GlobalSettings, ModelAdminBasics)
admin.site.register(models.Lab, LabAdmin)
admin.site.register(models.LabHead, ModelAdminBasics)
admin.site.register(models.LabProject)
admin.site.register(models.LabUserSettingsOverride, ModelAdminBasics)
admin.site.register(models.OrganizationUserSettingsOverride, ModelAdminBasics)
admin.site.register(models.Organization, OrganizationAdmin)
admin.site.register(models.Project, ModelAdminBasics)
admin.site.register(models.Sample, ModelAdminBasics)
admin.site.register(models.SampleLabProject)
admin.site.register(models.SampleTag, ModelAdminBasics)
admin.site.register(models.SettingsInitialGroupPermission, ModelAdminBasics)
admin.site.register(models.SiteMessage, ModelAdminBasics)
admin.site.register(models.Tag, ModelAdminBasics)
admin.site.register(models.Trio, ModelAdminBasics)
admin.site.register(models.UserSettingsOverride, UserSettingsOverrideAdmin)
admin.site.register(models.VCF, GuardedModelAdminBasics)
admin.site.register(models.VCFSourceSettings, ModelAdminBasics)
admin.site.register(models.VCFTag, ModelAdminBasics)
admin.site.register(models.VariantGridColumn, ModelAdminBasics)
