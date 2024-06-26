import re

from django.contrib import admin, messages
from django.contrib.admin.widgets import AdminTextInputWidget
from django.db.models import QuerySet
from unidecode import unidecode

from snpdb import models
from snpdb.admin_utils import ModelAdminBasics, GuardedModelAdminBasics, admin_list_column, \
    admin_action
from snpdb.liftover import liftover_alleles
from snpdb.models import Allele, VariantAllele, ClinVarKey, ClinVarKeyExcludePattern, UserSettingsOverride, \
    LabUserSettingsOverride, OrganizationUserSettingsOverride, UserPageAck, Organization, Lab, GlobalSettings, Variant, \
    AlleleLiftover
from snpdb.models.models_genome import GenomeBuild


@admin.register(AlleleLiftover)
class AlleleLiftoverAdmin(ModelAdminBasics):
    list_display = ("pk", "allele", "error_tidy_message", "status")
    list_filter = ("status", )
    search_fields = ('pk', 'allele__id')

    @admin_list_column("Error", order_field='error__message')
    def error_tidy_message(self, obj: AlleleLiftover):
        return obj.error_tidy()


@admin.register(Variant)
class VariantAdmin(ModelAdminBasics):
    pass


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
        genome_builds: list[GenomeBuild] = []
        for va in VariantAllele.objects.filter(allele=obj).order_by('genome_build'):
            genome_builds.append(va.genome_build)
        if genome_builds:
            return ", ".join(str(genome_build) for genome_build in genome_builds)
        return "-"

    @admin_action("Liftover")
    def liftover(self, request, queryset):
        liftover_alleles(allele_qs=queryset, user=request.user)
        self.message_user(request, message='Liftover queued', level=messages.INFO)


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


class SettingsOverrideAdmin(ModelAdminBasics):
    pass

    # Not sure hwo to just set timezone to a dropdown of available timezones, just need to be careful to enter an exact value
    # def get_form(self, request, obj=None, **kwargs):
    #
    #     return super().get_form(request, obj, widgets={
    #         'timezone': admin.widgets.AutocompleteSelect(choices=[(None, "")] + [(tz, tz) for tz in settings.AVAILABLE_TZS])
    #     }, **kwargs)


@admin.register(UserSettingsOverride)
class UserSettingsOverrideAdmin(ModelAdminBasics):
    list_display = ('user', 'default_lab', 'default_genome_build', 'email_weekly_updates', 'email_discordance_updates')
    list_filter = (DefaultBuildFilter,)
    search_fields = ('user__username', 'user__first_name', 'user__last_name')
    actions = ["export_as_csv"]

    def is_readonly_field(self, f):
        if f.name in ('default_genome_build', 'default_lab'):
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


@admin.register(GlobalSettings)
class GlobalSettingsAdmin(SettingsOverrideAdmin):
    pass


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
    list_display = ('id', 'name', 'last_full_run')
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

    @admin_action("ClinVar Export Prepare")
    def prepare_clinvar(self, request, queryset):
        from classification.models.clinvar_export_prepare import ClinvarExportPrepare
        report = ClinvarExportPrepare().update_export_records_for_keys(set(queryset.all()))
        if len(report) > 10:
            self.message_user(request, message="Showing first 10 messages", level=messages.INFO)
        for report_row in report[0:10]:
            self.message_user(request, message=report_row, level=messages.INFO)

    def get_form(self, request, obj=None, **kwargs):
        form = super().get_form(request, obj, widgets={
            'id': admin.widgets.AdminTextInputWidget(),
            'api_key': admin.widgets.AdminTextInputWidget(),
            'org_id': admin.widgets.AdminTextInputWidget(),
            'name': admin.widgets.AdminTextInputWidget()
        }, **kwargs)
        form.base_fields["assertion_method_lookup"].help_text = \
            'Preferred format is<br/>{"lookups":[<br/>&nbsp;{"match": "(regex1)", "citation": "acmg"},<br/>&nbsp;{"match": "(regex2)", "citation": {"db": "PubMed", "id": "PMID:123456"}}<br/>]}'
        return form


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
    search_fields = ('organization__name', 'name')

    fieldsets = (
        ('Basic', {'fields': ('name', 'group_name', 'organization')}),
        ('Position', {'fields': ('city', 'state', 'country', 'lat', 'long')}),
        ('Contact', {'fields': ('url', 'contact_name', 'contact_email', 'contact_phone')}),
        ('Notifications', {'fields': ('email', 'slack_webhook')}),
        ('Uploads', {'fields': ('upload_location', 'upload_automatic', 'upload_instructions')}),
        ('Submissions', {'fields': ('classification_config', 'external', 'clinvar_key', 'consolidates_variant_classifications')})
    )

    def is_readonly_field(self, f) -> bool:
        if f.name in ('clinvar_key', 'organization', 'state', 'country'):
            return False
        return super().is_readonly_field(f)

    def get_form(self, request, obj=None, **kwargs):

        return super().get_form(request, obj, widgets={
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
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'short_name': admin.widgets.AdminTextInputWidget(),
            'group_name': admin.widgets.AdminTextInputWidget()

        }, **kwargs)


@admin.register(models.LabHead)
class LabHeadAdmin(ModelAdminBasics):
    autocomplete_fields = (
        "lab",
        "user",
    )


@admin.register(models.UserAward)
class UserAwardAdmin(ModelAdminBasics):
    autocomplete_fields = (
        "user",
    )

    list_display = "user", "award_level", "award_text", "active"


@admin.register(models.Country)
class CountryAdmin(ModelAdminBasics):

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'short_name': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


@admin.register(models.State)
class StateAdmin(ModelAdminBasics):

    list_display = ('name', 'country')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'short_name': admin.widgets.AdminTextInputWidget()
        }, **kwargs)

    def is_readonly_field(self, f) -> bool:
        if f.name == "country":
            return False
        return super().is_readonly_field(f)


@admin.register(models.GenomicIntervalsCollection)
class GenomicIntervalsCollectionAdmin(ModelAdminBasics):
    search_fields = ('name', )


admin.site.register(models.CachedGeneratedFile, ModelAdminBasics)
admin.site.register(models.Cohort, ModelAdminBasics)
admin.site.register(models.CohortGenotypeCollection, ModelAdminBasics)
admin.site.register(models.CohortSample, ModelAdminBasics)
admin.site.register(models.CustomColumn, ModelAdminBasics)
admin.site.register(models.CustomColumnsCollection, ModelAdminBasics)
admin.site.register(models.GenomeBuild, ModelAdminBasics)
admin.site.register(models.LabProject)
admin.site.register(models.Manufacturer)
admin.site.register(models.Project, ModelAdminBasics)
admin.site.register(models.Sample, ModelAdminBasics)
admin.site.register(models.SampleLabProject)
admin.site.register(models.SampleTag, ModelAdminBasics)
admin.site.register(models.SettingsInitialGroupPermission, ModelAdminBasics)
admin.site.register(models.SiteMessage, ModelAdminBasics)
admin.site.register(models.Tag, ModelAdminBasics)
admin.site.register(models.Trio, ModelAdminBasics)
admin.site.register(models.VCF, GuardedModelAdminBasics)
admin.site.register(models.VCFSourceSettings, ModelAdminBasics)
admin.site.register(models.VCFTag, ModelAdminBasics)
admin.site.register(models.VariantGridColumn, ModelAdminBasics)
