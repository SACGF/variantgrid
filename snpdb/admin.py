import re
from typing import List

from django.contrib import admin, messages
from django.contrib.admin.widgets import AdminTextInputWidget
from django.db.models import QuerySet
from unidecode import unidecode

from snpdb import models
from snpdb.admin_utils import ModelAdminBasics, GuardedModelAdminBasics, admin_list_column, \
    admin_action
from snpdb.liftover import liftover_alleles
from snpdb.models import Allele, VariantAllele, ClinVarKey, ClinVarKeyExcludePattern, UserSettingsOverride, \
    LabUserSettingsOverride, OrganizationUserSettingsOverride, UserPageAck, Organization, Lab
from snpdb.models.models_genome import GenomeBuild


class AlleleClingenFilter(admin.SimpleListFilter):
    title = 'Clingen Filter'
    parameter_name = 'clin'
    default_value = None

    def lookups(self, request, model_admin):
        return [("present", "present"), ("missing", "missing")]

    def queryset(self, request, queryset):
        if value := self.value():
            if value == "present":
                return queryset.filter(clingen_allele__isnull=False)
            if value == "missing":
                return queryset.filter(clingen_allele__isnull=True)
        return queryset


@admin.register(Allele)
class AlleleAdmin(ModelAdminBasics):

    list_display = ('pk', 'clingen_allele', 'variants')
    list_filter = (AlleleClingenFilter,)
    search_fields = ('id', 'clingen_allele__id')

    @admin_list_column()
    def variants(self, obj: Allele):
        genome_builds: List[GenomeBuild] = list()
        for va in VariantAllele.objects.filter(allele=obj).order_by('genome_build'):
            genome_builds.append(va.genome_build)
        if genome_builds:
            return ", ".join(str(genome_build) for genome_build in genome_builds)
        return "-"

    @admin_action("Validate")
    def validate(self, request, queryset):
        allele: Allele
        for allele in queryset:
            allele.validate()

    @admin_action("Liftover")
    def liftover(self, request, queryset):
        liftover_alleles(allele_qs=queryset, user=request.user)
        self.message_user(request, message='Liftover queued', level=messages.INFO)

    @admin_action("ClinVar Export Prepare")
    def prepare_clinvar(self, request, queryset):
        from classification.models.clinvar_export_prepare import ClinvarAlleleExportPrepare
        allele: Allele
        for allele in queryset:
            export_prepare = ClinvarAlleleExportPrepare(allele=allele)
            report = export_prepare.update_export_records()
            for message in report:
                self.message_user(request, message=f"Allele ({allele}) - {message}", level=messages.INFO)


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


@admin.register(UserSettingsOverride)
class UserSettingsOverrideAdmin(ModelAdminBasics):
    list_display = ('user', 'default_lab', 'default_genome_build', 'email_weekly_updates', 'email_discordance_updates')
    list_filter = (DefaultBuildFilter,)
    search_fields = ('user__username', 'user__first_name', 'user__last_name')
    actions = ["export_as_csv"]

    def is_readonly_field(self, f):
        if f.name == 'default_genome_build':
            return False
        return super().is_readonly_field(f)


@admin.register(LabUserSettingsOverride)
class LabUserSettingsOverrideAdmin(ModelAdminBasics):
    list_display = ('lab', 'default_genome_build', 'email_weekly_updates', 'email_discordance_updates')
    actions = ["export_as_csv"]

    def is_readonly_field(self, f):
        if f.name == 'default_genome_build':
            return False
        return super().is_readonly_field(f)


@admin.register(OrganizationUserSettingsOverride)
class OrganizationUserSettingsOverrideAdmin(ModelAdminBasics):
    list_display = ('organization', 'default_genome_build', 'email_weekly_updates', 'email_discordance_updates')
    actions = ["export_as_csv"]

    def is_readonly_field(self, f):
        if f.name == 'default_genome_build':
            return False
        return super().is_readonly_field(f)


@admin.register(UserPageAck)
class UserPageAckAdmin(ModelAdminBasics):
    list_display = ('user', 'page_id')


class ClinVarKeyExcludePatternAdmin(admin.TabularInline):
    model = ClinVarKeyExcludePattern

    formfield_overrides = {
        models.TextField: {'widget': AdminTextInputWidget}
    }


@admin.register(ClinVarKey)
class ClinVarKeyAdmin(ModelAdminBasics):
    list_display = ('id', 'name', 'created', 'modified')
    inlines = (ClinVarKeyExcludePatternAdmin,)

    def run_ignores(self, request, queryset: QuerySet[ClinVarKey], apply: bool):
        from classification.models.clinvar_export_exclude_utils import ClinVarExcludePatternUtil
        for clinvar_key in queryset:
            ignore_util = ClinVarExcludePatternUtil(clinvar_key)
            test_results = ignore_util.run_all(apply=apply)
            if not apply:
                self.message_user(request, "Dry-run Only")
            if test_results:
                for key, ids in test_results.items():
                    message = f"{clinvar_key} {key} : {len(ids)} - e.g. classification ids {ids[0:5]}"
                    self.message_user(request, message)
            else:
                self.message_user(request, f"{clinvar_key} no ignores")

    @admin_action("Don't Share Flags - Dry Run")
    def ignore_dry_run(self, request, queryset: QuerySet[ClinVarKey]):
        self.run_ignores(request, queryset, False)

    @admin_action("Don't Share Flags - Apply")
    def ignore_apply(self, request, queryset: QuerySet[ClinVarKey]):
        self.run_ignores(request, queryset, True)

    def get_form(self, request, obj=None, **kwargs):

        return super(ClinVarKeyAdmin, self).get_form(request, obj, widgets={
            'id': admin.widgets.AdminTextInputWidget(),
            'api_key': admin.widgets.AdminTextInputWidget(),
            'org_id': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


def make_code_friendly(text: str) -> str:
    """
    convert accented characters to non-accented counterparts
    lower case, replace - and spaces with underscores
    remove anything that's then not a-z or underscore
    """
    text = unidecode(text) \
        .lower() \
        .replace('-', '_').replace(' ', '_')
    return re.sub(r'[^a-z0-9_]', '', text)


@admin.register(Lab)
class LabAdmin(ModelAdminBasics):
    list_per_page = 200
    list_display = ('name', 'group_name', 'organization', 'state', 'country',
                    'external', 'clinvar_key', 'upload_location', 'classification_config')

    fieldsets = (
        ('Basic', {'fields': ('name', 'group_name', 'organization')}),
        ('Position', {'fields': ('city', 'state', 'country', 'lat', 'long')}),
        ('Contact', {'fields': ('url', 'contact_name', 'contact_email', 'contact_phone')}),
        ('Notifications', {'fields': ('email', 'slack_webhook')}),
        ('Submissions', {'fields': ('classification_config', 'upload_location', 'upload_automatic', 'upload_instructions', 'external', 'clinvar_key')})
    )

    def is_readonly_field(self, f) -> bool:
        if f.name in ('clinvar_key', 'organization', 'state', 'country'):
            return False
        return super().is_readonly_field(f)

    def get_form(self, request, obj=None, **kwargs):

        return super(LabAdmin, self).get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'group_name': admin.widgets.AdminTextInputWidget(),
            'city': admin.widgets.AdminTextInputWidget(),
            'lat': admin.widgets.AdminTextInputWidget(),
            'long': admin.widgets.AdminTextInputWidget(),
            'url': admin.widgets.AdminURLFieldWidget(),
            'css_class': admin.widgets.AdminTextInputWidget(),
            'contact_name': admin.widgets.AdminTextInputWidget(),
            'contact_email': admin.widgets.AdminEmailInputWidget(),
            'contact_phone': admin.widgets.AdminTextInputWidget(),
            'upload_location': admin.widgets.AdminTextInputWidget(),
            'email': admin.widgets.AdminEmailInputWidget(),
            'slack_webhook': admin.widgets.AdminTextInputWidget()
        }, **kwargs)

    def fix_group_name(self, request, queryset):
        safety_reg = re.compile(r'^[a-z0-9_]*$')
        fixed = 0
        already_good = 0
        lab: Lab
        for lab in queryset:
            org_group_name = lab.organization.group_name
            if not lab.group_name or not safety_reg.match(lab.group_name) or not lab.group_name.startswith(org_group_name):
                lab.group_name = org_group_name + '/' + make_code_friendly(lab.name)
                lab.save()
                fixed = fixed + 1
            else:
                already_good = already_good + 1
        self.message_user(request, f"{fixed} updated, {already_good} no change required")

    fix_group_name.short_description = 'Fix group name'

    actions = [fix_group_name]


@admin.register(Organization)
class OrganizationAdmin(ModelAdminBasics):

    list_display = ('name', 'group_name', 'classification_config')

    fieldsets = (
        ('Basic', {'fields': ('name', 'short_name', 'group_name', 'active')}),
        ('Submissions', {'fields': ('classification_config', )})
    )

    def fix_group_name(self, request, queryset):
        org: Organization
        safety_reg = re.compile(r'^[a-z0-9_]*$')
        fixed = 0
        already_good = 0
        for org in queryset:
            if not org.group_name or not safety_reg.match(org.group_name):
                org.group_name = make_code_friendly(org.name)
                org.save()
                fixed = fixed + 1
            else:
                already_good = already_good + 1
        self.message_user(request, f"{fixed} updated, {already_good} no change required")

    fix_group_name.short_description = 'Fix group name'

    actions = [fix_group_name]

    def get_form(self, request, obj=None, **kwargs):
        return super(OrganizationAdmin, self).get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'short_name': admin.widgets.AdminTextInputWidget(),
            'group_name': admin.widgets.AdminTextInputWidget()

        }, **kwargs)


admin.site.register(models.CachedGeneratedFile, ModelAdminBasics)
admin.site.register(models.Cohort, ModelAdminBasics)
admin.site.register(models.CohortGenotypeCollection, ModelAdminBasics)
admin.site.register(models.CohortSample, ModelAdminBasics)
admin.site.register(models.Country)
admin.site.register(models.CustomColumn, ModelAdminBasics)
admin.site.register(models.CustomColumnsCollection, ModelAdminBasics)
admin.site.register(models.GenomeBuild, ModelAdminBasics)
admin.site.register(models.GenomicIntervalsCollection, ModelAdminBasics)
admin.site.register(models.GlobalSettings, ModelAdminBasics)
admin.site.register(models.LabHead, ModelAdminBasics)
admin.site.register(models.LabProject)
admin.site.register(models.Manufacturer)
admin.site.register(models.Project, ModelAdminBasics)
admin.site.register(models.Sample, ModelAdminBasics)
admin.site.register(models.SampleLabProject)
admin.site.register(models.SampleTag, ModelAdminBasics)
admin.site.register(models.SettingsInitialGroupPermission, ModelAdminBasics)
admin.site.register(models.SiteMessage, ModelAdminBasics)
admin.site.register(models.State)
admin.site.register(models.Tag, ModelAdminBasics)
admin.site.register(models.Trio, ModelAdminBasics)
admin.site.register(models.VCF, GuardedModelAdminBasics)
admin.site.register(models.VCFSourceSettings, ModelAdminBasics)
admin.site.register(models.VCFTag, ModelAdminBasics)
admin.site.register(models.VariantGridColumn, ModelAdminBasics)
