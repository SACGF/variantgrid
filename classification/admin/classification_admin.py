import json
from datetime import timedelta
from typing import Union, Optional

from django.contrib import admin, messages
from django.contrib.admin import RelatedFieldListFilter, BooleanFieldListFilter, DateFieldListFilter, TabularInline
from django.db.models import QuerySet
from django.forms import Widget
from django.utils import timezone
from django.utils.safestring import SafeString

from annotation.models.models import AnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.classification_import import reattempt_variant_matching, variant_matching_dry_run
from classification.enums import WithdrawReason
from classification.enums.classification_enums import EvidenceCategory, SpecialEKeys, SubmissionSource, ShareLevel
from classification.models import EvidenceKey, EvidenceKeyMap, DiscordanceReport, DiscordanceReportClassification, \
    ClinicalContext, ClassificationReportTemplate, ClassificationModification, \
    UploadedClassificationsUnmapped, ImportedAlleleInfo, ClassificationImport, ImportedAlleleInfoStatus, \
    classification_flag_types, DiscordanceReportTriage, ensure_discordance_report_triages_bulk, \
    DiscordanceReportTriageStatus, ClassificationGrouping, ClassificationGroupingEntry, \
    AlleleOriginGrouping, AlleleGrouping, ClassificationGroupingSearchTerm
from classification.models.classification import Classification
from classification.models.classification_import_run import ClassificationImportRun, ClassificationImportRunStatus
from classification.models.classification_variant_info_models import ResolvedVariantInfo, ImportedAlleleInfoValidation
from classification.models.clinical_context_models import ClinicalContextRecalcTrigger, DiscordanceNotification
from classification.models.clinical_context_utils import update_clinical_contexts
from classification.models.discordance_lab_summaries import DiscordanceLabSummary
from classification.models.discordance_models_utils import DiscordanceReportRowDataTriagesRowData
from classification.signals import send_prepared_discordance_notifications
from classification.tasks.classification_import_map_and_insert_task import ClassificationImportMapInsertTask
from library.cache import timed_cache
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from library.utils import ExportRow, export_column, ExportDataType, first
from ontology.models import OntologyTerm, AncestorCalculator
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column, AllValuesChoicesFieldListFilter, \
    admin_model_action, get_admin_url
from snpdb.lab_picker import LabPickerData
from snpdb.models import GenomeBuild, Lab


class VariantMatchedFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'variant match status'
    parameter_name = 'matched'
    default_value = None

    def lookups(self, request, model_admin):
        return [('True', 'Matched'), ('False', 'Not-Matched')]

    def queryset(self, request, queryset):
        if self.value() == 'True':
            return queryset.filter(variant__isnull=False)
        if self.value() == 'False':
            return queryset.filter(variant__isnull=True)
        return queryset


class ClassificationShareLevelFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'share level'
    parameter_name = 'share_level'
    default_value = None

    def lookups(self, request, model_admin):
        return [(sl.key, sl.label) for sl in ShareLevel.ALL_LEVELS]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(share_level=self.value())
        return queryset


class ClinicalContextFilter(admin.SimpleListFilter):
    title = 'clinical context'
    parameter_name = 'clinical_context'
    default_value = 'default'

    def lookups(self, request, model_admin):
        return [('missing', 'Missing'), ('default', 'Default'), ('custom', 'Custom')]

    def queryset(self, request, queryset):
        lookup = self.value()
        if lookup == 'missing':
            return queryset.filter(clinical_context__isnull=True)
        if lookup == 'default':
            return queryset.filter(clinical_context__name='default')
        if lookup == 'custom':
            return queryset.filter(clinical_context__isnull=False).exclude(clinical_context__name='default')
        return queryset


class ClassificationImportedGenomeBuildFilter(admin.SimpleListFilter):
    title = 'imported genome build'
    parameter_name = 'genome_build'
    default_value = None

    def lookups(self, request, model_admin):
        builds = GenomeBuild.objects.all()
        return [(build.pk, build.name) for build in builds]

    def queryset(self, request, queryset):
        if value := self.value():
            substring = GenomeBuild.objects.get(pk=value).name
            return queryset.filter(evidence__genome_build__value__istartswith=substring)
        return queryset


class ClassificationLatestImportFilter(admin.SimpleListFilter):
    title = "import (choose a lab first)"
    parameter_name = "latest_import"
    default_value = None

    def lookups(self, request, model_admin):
        return [("latest", "latest"), ("not-latest", "not-latest")]

    def queryset(self, request, queryset: QuerySet[Classification]):
        if value := self.value():
            most_up_to_date_run: Optional[ClassificationImportRun] = None
            if most_up_to_date_class_run := queryset.filter(last_import_run__isnull=False).order_by(
                    '-last_import_run__created').first():
                most_up_to_date_run = most_up_to_date_class_run.last_import_run
            if most_up_to_date_run:
                if value == "latest":
                    return queryset.filter(last_import_run=most_up_to_date_run)
                if value == "not-latest":
                    return queryset.exclude(last_import_run=most_up_to_date_run)
        return queryset


class ClassificationModificationAdmin(admin.TabularInline):
    model = ClassificationModification

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


class ClassificationImportRunCountFilter(admin.SimpleListFilter):
    title = "Run Type"
    parameter_name = "upload_type"
    default_value = None

    def lookups(self, request, model_admin):
        return [
            ("import", "Proper Import"),
            ("test", "Import Test Tool"),
            ("rematch", "Allele Re-Match"),
        ]

    def queryset(self, request, queryset: QuerySet[ClassificationImportRun]):
        if self.value() == "test":
            return queryset.filter(identifier__endswith="#classification_import_tool")
        elif self.value() == "rematch":
            return queryset.filter(identifier="admin-variant-rematch")
        elif self.value() == "import":
            return queryset.exclude(identifier__endswith="#classification_import_tool").exclude(identifier="admin-variant-rematch")


@admin.register(ClassificationImportRun)
class ClassificationImportRunAdmin(ModelAdminBasics):
    # change_list_template = 'classification/admin/change_list.html'
    search_fields = ('identifier', 'from_file__filename',)
    list_display = ['id', 'identifier', 'row_count', 'status', 'from_file', 'created_detailed', 'modified_detailed']
    list_filter = ('status', ClassificationImportRunCountFilter)

    @admin_model_action(url_slug="create_dummy/", short_description="Create Dummy Import", icon="fa-solid fa-plus")
    def create_dummy(self, request):
        ClassificationImportRun(identifier="Dummy Import").save()
        self.message_user(request, "Created a dummy import")

    @admin_model_action(url_slug="mark_unfinished/", short_description="Mark OnGoing Imports as Unfinished",
                        icon="fa-regular fa-trash-can")
    def close_all_open(self, request):
        ongoing_imports = ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING)
        if ongoing_imports:
            for cir in ongoing_imports:
                cir.status = ClassificationImportRunStatus.UNFINISHED
                cir.save()
                self.message_user(request, message=f'Changed import from ongoing to unfinished for {cir.identifier}',
                                  level=messages.INFO)
        else:
            self.message_user(request, message='No OnGoing Imports to mark as unfinished',
                              level=messages.WARNING)

    def is_readonly_field(self, f) -> bool:
        return True

    @admin_list_column(short_description="Created", order_field="created")
    def created_detailed(self, obj: ClassificationImportRun):
        return self.format_datetime(obj.created)

    @admin_list_column(short_description="Modified", order_field="modified")
    def modified_detailed(self, obj: ClassificationImportRun):
        return self.format_datetime(obj.modified)

    @admin_action(short_description="Mark Unfinished")
    def mark_unfinished(self, request, queryset: QuerySet[ClassificationImportRun]):
        for cir in queryset:
            if cir.status == ClassificationImportRunStatus.ONGOING:
                cir.status = ClassificationImportRunStatus.UNFINISHED
                cir.save()
                self.message_user(request, message='Changed import from ongoing to unfinished',
                                  level=messages.INFO)
            else:
                self.message_user(request, message='Can only mark unfinished ongoing records',
                                  level=messages.WARNING)


@admin.register(Classification)
class ClassificationAdmin(ModelAdminBasics):
    list_display = [
        'id',
        'lab',
        'lab_record_id',
        'last_import_run',
        'last_source_id',
        'share_level',
        'clinical_significance',
        'allele_fallback',
        'grch37_c_hgvs',
        'grch38_c_hgvs',
        'imported_genome_build',
        'imported_c_hgvs',
        'withdrawn',
        'user',
        'created_detailed',
        'modified_detailed']
    list_filter = (
        ('lab__organization', RelatedFieldListFilter),
        ('lab', RelatedFieldListFilter),
        ClassificationLatestImportFilter,
        ('withdrawn', BooleanFieldListFilter),
        ClassificationShareLevelFilter,
        VariantMatchedFilter,
        ClinicalContextFilter,
        ClassificationImportedGenomeBuildFilter,
        'allele_origin_bucket'
        # ('user', RelatedFieldListFilter),
    )
    search_fields = ('id', 'lab_record_id')
    list_per_page = 50
    inlines = (ClassificationModificationAdmin,)
    list_select_related = ('lab', 'user', 'allele')

    @admin_list_column(short_description="c.hgvs (37)", order_field="allele_info__grch37__c_hgvs")
    def grch37_c_hgvs(self, obj: Classification):
        try:
            return obj.allele_info.grch37.c_hgvs
        except AttributeError:
            return ""

    @admin_list_column(short_description="c.hgvs (38)", order_field="allele_info__grch37__c_hgvs")
    def grch38_c_hgvs(self, obj: Classification):
        try:
            return obj.allele_info.grch38.c_hgvs
        except AttributeError:
            return ""

    @admin_list_column(short_description="Created", order_field="created")
    def created_detailed(self, obj: Classification):
        return self.format_datetime(obj.created)

    @admin_list_column(short_description="Modified", order_field="modified")
    def modified_detailed(self, obj: Classification):
        return self.format_datetime(obj.modified)

    @admin_list_column(short_description="Allele", order_field="allele")
    def allele_fallback(self, obj: Classification):
        if allele := obj.allele:
            return str(allele)
        if variant := obj.variant:
            return str(variant)
        return "-"

    def has_add_permission(self, request):
        return False

    @admin_action("Data: Populate w Annotations")
    def populate_base_variant_data(self, request, queryset: QuerySet[Classification]):
        for vc in queryset:
            refseq_transcript_id = vc.get(SpecialEKeys.REFSEQ_TRANSCRIPT_ID)
            ensembl_transcript_id = vc.get(SpecialEKeys.ENSEMBL_TRANSCRIPT_ID)
            if not vc.variant:
                self.message_user(request, "(%i) No variant for classification" % vc.id)
            else:
                genome_build = vc.get_genome_build()
                annotation_version = AnnotationVersion.latest(genome_build)
                data = get_evidence_fields_for_variant(genome_build, vc.variant, refseq_transcript_id,
                                                       ensembl_transcript_id,
                                                       evidence_keys_list=[], annotation_version=annotation_version)
                patch = {}
                publish = False
                for key in [SpecialEKeys.GENOME_BUILD,
                            SpecialEKeys.G_HGVS,
                            SpecialEKeys.C_HGVS,
                            SpecialEKeys.P_HGVS,
                            SpecialEKeys.VARIANT_COORDINATE,
                            SpecialEKeys.CLINGEN_ALLELE_ID]:
                    if key in data:
                        patch[key] = data[key]

                response = vc.patch_value(
                    patch=patch,
                    source=SubmissionSource.VARIANT_GRID,
                    save=True,
                    user=admin_bot(),
                    leave_existing_values=True
                )
                actual_patch = {}
                if response.get('modified'):
                    for key in response.get('modified'):
                        actual_patch[key] = patch[key]
                    vc.publish_latest(user=admin_bot())
                    publish = True

                if publish:
                    self.message_user(request, "(%i) Patched = %s" % (vc.id, json.dumps(actual_patch)))
                else:
                    self.message_user(request, "(%i) No changes with = %s" % (vc.id, json.dumps(patch)))

    @admin_action("Fixes: Condition Text Sync")
    def condition_text_sync(self, request, queryset: QuerySet[Classification]):
        from classification.models import ConditionTextMatch

        conditions_newly_set = 0
        for c in queryset:
            has_condition_already = bool(c.condition_resolution)
            ConditionTextMatch.sync_condition_text_classification(c.last_published_version)
            c.refresh_from_db(fields=["condition_resolution"])
            has_condition_now = bool(c.condition_resolution)
            if has_condition_now and not has_condition_already:
                conditions_newly_set += 1

        self.message_user(request, f"{conditions_newly_set} records have conditions now when they didn't previously")

    @admin_action("Fixes: Fix Permissions")
    def fix_permissions(self, request, queryset: QuerySet[Classification]):
        for c in queryset:
            c.fix_permissions(fix_modifications=True)

    @admin_action("Fixes: Clinical Context Germline/Somatic")
    def fix_clinical_context(self, request, queryset: QuerySet[Classification]):
        update_clinical_contexts(list(queryset.all()))

    @admin_action("Fixes: Revalidate")
    def revalidate(self, request, queryset):
        for vc in queryset:
            vc.revalidate(request.user)
        self.message_user(request, str(queryset.count()) + " records revalidated")

    @admin_action("Matching: Re-Match Variant")
    def admin_reattempt_variant_matching(self, request, queryset: QuerySet[Classification]):
        for classification in queryset:
            _, created = classification.ensure_allele_info_with_created(force_allele_info_update_check=True)
            if created:
                classification.save()

        allele_info_ids = queryset.values_list('allele_info', flat=True)
        allele_info_qs = ImportedAlleleInfo.objects.filter(id__in=allele_info_ids)

        reattempt_variant_matching(request.user, allele_info_qs, False)

    @admin_action("Matching: Re-Calculate c.hgvs, transcripts")
    def recalculate_cached_chgvs(self, request, queryset: QuerySet[Classification]):
        for vc in queryset:
            vc.update_cached_c_hgvs()
            vc.save()

    @admin_action("Matching: Re-apply AlleleInfo to Classification")
    def re_apply_allele_info(self, request, queryset: QuerySet[Classification]):
        for classification in queryset:
            classification.attempt_set_variant_info_from_pre_existing_imported_allele_info()
            classification.save()

    def publish_share_level(self, request, queryset: QuerySet[Classification], share_level: ShareLevel):
        already_published = 0
        in_error = 0
        published = 0
        for vc in queryset:
            if vc.has_errors():
                in_error += 1
            elif vc.share_level_enum >= share_level and not vc.has_outstanding_changes():
                already_published += 1
            else:
                try:
                    vc.publish_latest(user=request.user, share_level=share_level)
                    published += 1
                except ValueError as ve:
                    self.message_user(request, message='Error when publishing - ' + str(ve),
                                      level=messages.ERROR)
            if not in_error and not vc.has_outstanding_changes():
                vc.flag_collection.ensure_resolution(
                    flag_type=classification_flag_types.classification_outstanding_edits,
                    resolution='closed'
                )

        if already_published:
            self.message_user(request, message=f"({already_published}) records had been previously published",
                              level=messages.INFO)
        if in_error:
            self.message_user(request, message=f"({in_error}) records can't be published due to validation errors",
                              level=messages.ERROR)
        if published:
            self.message_user(request, message=f"({published}) records have been freshly published",
                              level=messages.INFO)

    @admin_action("Publish: Organisation")
    def publish_org(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.INSTITUTION)

    @admin_action("Publish: Logged in Users")
    def publish_logged_in_users(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.ALL_USERS)

    @admin_action("State: Make Mutable")
    def make_mutable(self, request, queryset: QuerySet[Classification]):
        for vc in queryset:
            vc.patch_value(patch={}, user=request.user, source=SubmissionSource.VARIANT_GRID, remove_api_immutable=True,
                           save=True)

    def set_withdraw(self, request, queryset: QuerySet[Classification], withdraw: bool,
                     reason: WithdrawReason = WithdrawReason.OTHER) -> int:
        count = 0
        for vc in queryset:
            try:
                actioned = vc.set_withdrawn(user=request.user, withdraw=withdraw, reason=reason)
                if actioned:
                    count += 1
            except BaseException:
                pass
        return count

    @admin_action("State: Withdraw (Other)")
    def withdraw_true_other(self, request, queryset: QuerySet[Classification]):
        count = self.set_withdraw(request, queryset, True)
        self.message_user(request, f"{count} records now newly set to withdrawn")

    @admin_action("State: Withdraw (Duplicate)")
    def withdraw_true_duplicate(self, request, queryset: QuerySet[Classification]):
        count = self.set_withdraw(request, queryset, True, reason=WithdrawReason.DUPLICATE)
        self.message_user(request, f"{count} records now newly set to withdrawn as duplicate")

    @admin_action("State: Un-Withdraw")
    def withdraw_false(self, request, queryset):
        count = self.set_withdraw(request, queryset, False)
        self.message_user(request, f"{count} records now newly set to un-withdrawn")

    """
    @admin_action("Fix allele Freq History")
    def fix_allele_freq_history(self, request, queryset):
        results = EvidenceKeyToUnit(key_names=["allele_frequency"]).migrate(queryset, dry_run=False)
        for result in results:
            self.message_user(request, result)
    """

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'lab_record_id': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


class ImportedAlleleInfoInline(admin.TabularInline):
    model = ImportedAlleleInfo
    fields = ['id', 'imported_c_hgvs', 'imported_genome_build_patch_version', 'status']

    def has_add_permission(self, request, obj):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


@admin.register(ClassificationImport)
class ClassificationImportAdmin(ModelAdminBasics):
    inlines = (ImportedAlleleInfoInline,)
    list_display = ('pk', 'created', 'user', 'genome_build', 'outstanding_classifications')

    @admin_list_column("Outstanding Allele Infos")
    def outstanding_classifications(self, obj: ClassificationImport):
        return ImportedAlleleInfo.objects.filter(classification_import=obj).count()

    @admin_action("Count Variants")
    def count_variants(self, request, queryset: QuerySet[ClassificationImport]):
        for obj in queryset:
            self.message_user(request, f'Classification Import ({obj.pk}) has {obj.get_variants_qs().count()} variants')

    @admin_action("Relink Variants")
    def relink_variants(self, request, queryset: QuerySet[ClassificationImport]):
        for obj in queryset:
            updated_ai_count = ImportedAlleleInfo.relink_variants(obj)
            self.message_user(request, f'Classification Import ({obj.pk}) relinked {updated_ai_count} variants')


@admin.register(ClinicalContext)
class ClinicalContextAdmin(ModelAdminBasics):
    list_display = (
    'id', 'allele', 'name', 'status', 'allele_origin_bucket', 'modified', 'pending_cause', 'pending_status')
    search_fields = ('id', 'allele__pk', 'name')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget()
        })

    def has_add_permission(self, request):
        return False

    @admin_action("Recalculate Status")
    def recalculate(self, request, queryset):
        for dc in queryset:
            dc.recalc_and_save(cause='Admin recalculation',
                               cause_code=ClinicalContextRecalcTrigger.ADMIN)  # cause of None should change to Unknown, which is accurate if this was required
        self.message_user(request, 'Recalculated %i statuses' % queryset.count())


class EvidenceKeySectionFilter(admin.SimpleListFilter):
    title = 'Section Filter'
    parameter_name = 'section'
    default_value = None

    def lookups(self, request, model_admin):
        return EvidenceCategory.CHOICES

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(evidence_category=self.value()).order_by('order', 'key')
        return queryset


class MaxShareLevelFilter(admin.SimpleListFilter):
    title = 'max share level'
    parameter_name = 'max_share_level'
    default_value = None

    def lookups(self, request, model_admin):
        return ShareLevel.choices()

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(max_share_level=self.value()).order_by('order', 'key')
        return queryset


@admin.register(EvidenceKey)
class EvidenceKeyAdmin(ModelAdminBasics):
    list_per_page = 500
    list_filter = (EvidenceKeySectionFilter, MaxShareLevelFilter)
    ordering = ('key',)
    list_display = ('key', 'label', 'value_type', 'max_share_level', 'evidence_category', 'order')
    search_fields = ('key', 'label')

    fieldsets = (
        ('Basic', {'fields': ('key', 'label', 'sub_label')}),
        ('Position', {'fields': ('hide', 'evidence_category', 'order')}),
        ('Type', {'fields': ('mandatory', 'value_type', 'options', 'allow_custom_values', 'default_crit_evaluation',
                             'crit_allows_override_strengths', 'crit_uses_points')}),
        ('Overrides', {'fields': ('namespace_overrides',)}),
        ('Help', {'fields': ('description', 'examples', 'see')}),
        ('Admin', {'fields': ('max_share_level', 'copy_consensus', 'variantgrid_column', 'immutable')})
    )

    @admin_action("Validate key")
    def validate(self, request, queryset):
        good_count = 0
        key: EvidenceKey
        for key in queryset:
            try:
                key.validate()
                good_count += 1
            except ValueError as ve:
                self.message_user(request, str(ve), 'error')

        self.message_user(request, str(good_count) + " keys passed validation")

    def widget_overrides(self) -> dict[str, Widget]:
        return {
            'key': admin.widgets.AdminTextInputWidget(),
            'label': admin.widgets.AdminTextInputWidget(),
            'sub_label': admin.widgets.AdminTextInputWidget(),
            'see': admin.widgets.AdminURLFieldWidget()
        }


@admin.register(ClassificationReportTemplate)
class ClassificationReportTemplateAdmin(admin.ModelAdmin):
    list_display = ('name', 'modified')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'template': admin.widgets.AdminTextareaWidget()
        }, **kwargs)


class DiscordanceReportAdminLabFilter(admin.SimpleListFilter):
    title = "Labs Involved"
    parameter_name = "labs_involved"

    def lookups(self, request, model_admin):
        return [("multilabs", "MultipleLabs")]

    def queryset(self, request, queryset):
        if self.value() == "multilabs":
            valid_ids = set()
            for obj in queryset:
                labs = set()
                for drc in obj.discordancereportclassification_set.all():
                    if c := drc.classification_original.classification:
                        labs.add(c.lab)
                if len(labs) > 1:
                    valid_ids.add(obj.id)

            queryset = DiscordanceReport.objects.filter(pk__in=valid_ids)
        return queryset


class DiscordanceReportClassificationAdmin(admin.TabularInline):
    model = DiscordanceReportClassification
    readonly_fields = ["classification_original", "clinical_context_effective", "clinical_context_final",
                       "withdrawn_final"]
    fields = ["classification_original", "clinical_context_effective", "clinical_context_final", "withdrawn_final"]

    def clinical_context_effective(self, drc: DiscordanceReportClassification):
        return drc.clinical_context_effective

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


class DiscordanceReportAdminExport(ExportRow):

    @staticmethod
    @timed_cache(ttl=60)
    def cs_to_index():
        return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("vg")

    @staticmethod
    def _less_more_certain(summary: DiscordanceLabSummary):
        cs_to_index = DiscordanceReportAdminExport.cs_to_index()
        from_value = int(cs_to_index.get(summary.clinical_significance_from, "0"))
        to_value = int(cs_to_index.get(summary.clinical_significance_to, "0"))

        if summary.clinical_significance_to == 'withdrawn':
            return "withdrawn"
        elif from_value == to_value:
            return "same"
        elif from_value == 0 or to_value == 0:
            return "?"
        else:
            if abs(to_value - 3) > abs(from_value - 3):
                return "more"
            else:
                return "less"

    @staticmethod
    def _up_down_for(summary: DiscordanceLabSummary):
        cs_to_index = DiscordanceReportAdminExport.cs_to_index()
        from_value = int(cs_to_index.get(summary.clinical_significance_from, "0"))
        to_value = int(cs_to_index.get(summary.clinical_significance_to, "0"))

        if summary.clinical_significance_to == 'withdrawn':
            return "withdrawn"
        elif from_value == to_value:
            return "same"
        elif from_value == 0 or to_value == 0:
            return "?"
        else:
            if to_value > from_value:
                return "upgrade"
            else:
                return "downgrade"

    def __init__(self, discordance_report: DiscordanceReport, perspective: LabPickerData):
        self.discordance_report = discordance_report
        self.summaries = DiscordanceLabSummary.for_discordance_report(discordance_report, perspective)

    @export_column("Report number")
    def _report_number(self):
        return self.discordance_report.pk

    @export_column("URL")
    def _url(self):
        return get_url_from_view_path(self.discordance_report.get_absolute_url())

    @export_column("# labs")
    def _labs_involved(self):
        return len({summary.lab for summary in self.summaries})

    @export_column("# classifications")
    def _classification_count(self):
        return sum((summary.count for summary in self.summaries))

    @export_column("Status")
    def _status(self):
        return self.discordance_report.resolution_text or "Unknown"

    @export_column("Date Opened", data_type=ExportDataType.date)
    def _date_opened(self):
        return self.discordance_report.report_started_date

    @export_column("Date Closed", data_type=ExportDataType.date)
    def _date_closed(self):
        return self.discordance_report.report_completed_date

    @export_column("Days in Discordance")
    def _days_in_discordance(self):
        dr = self.discordance_report
        started = dr.report_started_date
        delta: timedelta
        if closed := dr.report_completed_date:
            delta = closed - started
        else:
            delta = timezone.now() - started

        days_delta_float = delta.seconds / 144.0
        if days_delta_float < 7:
            return f"{days_delta_float:0.1f}"
        else:
            return int(days_delta_float)

    @export_column("Gene Symbol")
    def _gene_symbol(self):
        all_chgvs = ImportedAlleleInfo.all_chgvs(self.discordance_report.clinical_context.allele)
        return "\n".join(sorted({chgvs.gene_symbol for chgvs in all_chgvs}))

    @export_column("c.HGVS (38)")
    def _variant(self):
        all_chgvs = ImportedAlleleInfo.all_chgvs(self.discordance_report.clinical_context.allele)
        c38s = sorted([str(chgvs) for chgvs in all_chgvs if chgvs.genome_build == GenomeBuild.grch38()])
        if c38s:
            return "\n".join(c38s)

    @export_column("Admin Notes")
    def _admin_notes(self):
        return self.discordance_report.admin_note

    @export_column("Umbrella Condition")
    def _umbrella_condition(self):
        try:
            umbrellas = ""
            conditions: set[OntologyTerm] = set()
            unresolved_conditions: set[str] = set()
            for summary in self.summaries:
                for drc in summary.drcs:
                    has_resolved_conditions = False
                    if condition_obj := drc.classification_original.classification.condition_resolution_obj:
                        if mondo := condition_obj.as_mondo_if_possible():
                            conditions.update(mondo.terms)
                            has_resolved_conditions = True
                    if not has_resolved_conditions:
                        if raw_condition := drc.classification_effective.get(SpecialEKeys.CONDITION):
                            unresolved_conditions.add(raw_condition)

            if conditions:
                ancestor = AncestorCalculator.common_ancestor(conditions)
                umbrellas = str(ancestor)

            if unresolved_conditions:
                unresolved_str = ", ".join(sorted(unresolved_conditions))
                if umbrellas:
                    umbrellas += f", **({unresolved_str})"
                else:
                    umbrellas = unresolved_str

            return umbrellas.strip()
        except ValueError:
            return "No common ancestor could be found"

    @export_column("Clinically Significant")
    def _clinically_significant(self):
        for summary in self.summaries:
            if "P" in summary.clinical_significance_from:
                return True
        return False

    @export_column("Overall Upgrade/Downgrade")
    def _overall_upgrade_downgrade(self):
        ups_and_downs = {DiscordanceReportAdminExport._up_down_for(summary) for summary in self.summaries}
        ups_and_downs.discard("same")
        if not ups_and_downs:
            return "same"
        elif len(ups_and_downs) == 1:
            return first(ups_and_downs)
        if "upgrade" in ups_and_downs and "downgrade" in ups_and_downs:
            return "mixed"
        elif "upgrade" in ups_and_downs:
            return "upgrade"
        elif "downgrade" in ups_and_downs:
            return "downgrade"
        else:
            return first(ups_and_downs)

    @export_column("Overall Certainty")
    def _overall_certainty(self):
        certainties = {DiscordanceReportAdminExport._less_more_certain(summary) for summary in self.summaries}
        certainties.discard("same")
        if not certainties:
            return "same"
        elif len(certainties) == 1:
            return first(certainties)
        if "more" in certainties and "less" in certainties:
            return "mixed"
        elif "more" in certainties:
            return "more"
        elif "less" in certainties:
            return "less"
        else:
            return first(certainties)

    @export_column("Labs")
    def _labs(self):
        return "\n".join(str(summary.lab) for summary in self.summaries)

    @export_column("Conditions")
    def _conditions(self):
        condition_rows = []
        for summary in self.summaries:
            conditions = {drc.classification_effective.condition_text or "" for drc in summary.drcs}
            condition_rows.append(", ".join(sorted(conditions)))
        return "\n".join(condition_rows)

    @export_column("Classification (Original)")
    def _cs_original(self):
        return "\n".join(str(summary.clinical_significance_from) for summary in self.summaries)

    @export_column("Classification (Current)")
    def _cs_current(self):
        return "\n".join(str(summary.clinical_significance_to) for summary in self.summaries)

    @export_column("Upgrade/Downgrade")
    def _upgrade_downgrade(self):
        return "\n".join((DiscordanceReportAdminExport._up_down_for(summary) for summary in self.summaries))

    @export_column("Certainty")
    def _certainty(self):
        return "\n".join((DiscordanceReportAdminExport._less_more_certain(summary) for summary in self.summaries))

    @export_column("Triage", sub_data=DiscordanceReportRowDataTriagesRowData)
    def _triages(self):
        return DiscordanceReportRowDataTriagesRowData(discordance_report=self.discordance_report,
                                                      perspective=LabPickerData.for_admin())


@admin.register(DiscordanceReport)
class DiscordanceReportAdmin(ModelAdminBasics):
    list_display = ["pk", "report_started_date", "c_hgvs", "days_open", "classification_count", "clinical_sigs", "labs",
                    "anotes"]
    list_select_related = ('clinical_context', 'clinical_context__allele')
    list_filter = [DiscordanceReportAdminLabFilter]
    inlines = (DiscordanceReportClassificationAdmin,)

    # @admin_list_column("Allele", order_field="clinical_context__allele__pk")
    # def allele(self, obj: DiscordanceReport) -> str:
    #     cc = obj.clinical_context
    #     return str(cc.allele)

    @admin_list_column("Admin Notes")
    def anotes(self, obj: DiscordanceReport):
        # make this an admin list column, so it crops the characters
        return obj.admin_note

    @admin_list_column("c.HGVS")
    def c_hgvs(self, obj: DiscordanceReport):
        for record in obj.discordancereportclassification_set.select_related(
                'classification_original__classification__allele_info'):
            if allele_info := record.classification_original.classification.allele_info:
                if c_hgvs := allele_info.grch38:
                    return c_hgvs
        # if can't find any 38 c.HGVSs, fallback onto Allele
        cc = obj.clinical_context
        return str(cc.allele.clingen_allele)

    @admin_list_column("Classifications")
    def clinical_sigs(self, obj: DiscordanceReport) -> str:
        clinical_sigs = set()
        for dr in DiscordanceReportClassification.objects.filter(report=obj):
            clinical_sigs.add(dr.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

        e_key_cs = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        sorter = e_key_cs.classification_sorter_value
        sorted_list = sorted(clinical_sigs, key=lambda x: (sorter(x), x))
        return sorted_list
        # pretty = ", ".join([e_key_cs.pretty_value(cs) for cs in sorted_list])
        # return pretty

    @admin_list_column()
    def days_open(self, obj: DiscordanceReport) -> Union[int, str]:
        closed = obj.report_completed_date
        if not closed:
            delta = timezone.now() - obj.report_started_date
            return f"{delta.days}+"
        delta = closed - obj.report_started_date
        return delta.days

    @admin_list_column()
    def classification_count(self, obj: DiscordanceReport):
        return obj.discordancereportclassification_set.count()

    @admin_list_column()
    def labs(self, obj: DiscordanceReport) -> str:
        labs: set[Lab] = set()
        drc: DiscordanceReportClassification
        for drc in obj.discordancereportclassification_set.all():
            if c := drc.classification_original.classification:
                labs.add(c.lab)
        labs_sorted = sorted(labs, key=lambda x: x.name)
        return ", ".join([lab.name for lab in labs_sorted])

    def has_add_permission(self, request):
        return False

    @admin_action("Re-calculate latest")
    def re_calculate(self, request, queryset):
        ds: DiscordanceReport
        for ds in queryset:
            ds.clinical_context.recalc_and_save(cause="Admin recalculation",
                                                cause_code=ClinicalContextRecalcTrigger.ADMIN)

    @admin_action("Export Admin Report CSV")
    def export_admin_report(self, request, queryset: QuerySet[DiscordanceReport]):
        perspective = LabPickerData.for_user(request.user)
        return DiscordanceReportAdminExport.streaming(request, (DiscordanceReportAdminExport(dr, perspective) for dr in
                                                                queryset), filename="discordance_admin_report")


class TriageStatusFilter(admin.SimpleListFilter):
    title = "Status"
    parameter_name = "status"

    def lookups(self, request, model_admin):
        return [("not_pending", "Not Pending")]

    def queryset(self, request, queryset: QuerySet[DiscordanceReportTriage]):
        if self.value() == "not_pending":
            queryset = queryset.exclude(triage_status=DiscordanceReportTriageStatus.PENDING)
        return queryset


@admin.register(DiscordanceReportTriage)
class DiscordanceReportTriageAdmin(ModelAdminBasics):
    list_display = ("pk", "lab", "triage_status_extra", "discordance_report_extra", "user", "modified")
    list_filter = (TriageStatusFilter, "closed")

    def is_readonly_field(self, f) -> bool:
        if f.name == "user":
            return False
        return super().is_readonly_field(f)

    @admin_model_action(url_slug="ensure_bulk/", short_description="Ensure Bulk", icon="fa-solid fa-dolly")
    def ensure_bulk(self, request):
        ensure_discordance_report_triages_bulk()

    @admin_list_column("Discordance Report", order_field="discordance_report", limit=100)
    def discordance_report_extra(self, obj: DiscordanceReportTriage):
        dr = obj.discordance_report
        total_triages = dr.discordancereporttriage_set.count()
        completed_triages = dr.discordancereporttriage_set.exclude(
            triage_status=DiscordanceReportTriageStatus.PENDING).count()
        url = get_admin_url(dr)
        return SafeString(
            f'(<a href="{url}">DR_{dr.pk}</a>) {dr.get_resolution_display() or "Discordant"} - <span style="font-size:0.6rem">triages complete {completed_triages} <i>of</i> {total_triages}</span>'
        )

    @admin_list_column("Triage Status", order_field="triage_status", limit=100)
    def triage_status_extra(self, obj: DiscordanceReportTriage):
        report = obj.get_triage_status_display()
        if obj.closed:
            report += " CLOSED"
        return report


@admin.register(DiscordanceNotification)
class DiscordanceNotificationAdmin(ModelAdminBasics):
    list_display = ("lab", "discordance_report", "notification_sent_date")
    list_filter = (('lab', RelatedFieldListFilter), ('notification_sent_date', DateFieldListFilter))

    def has_add_permission(self, request):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False

    @admin_action("Send Notification")
    def send_lab_discordance_notification(self, request, queryset):
        send_prepared_discordance_notifications(queryset)


@admin.register(UploadedClassificationsUnmapped)
class UploadedClassificationsUnmappedAdmin(ModelAdminBasics):
    list_display = ("pk", "lab", "created", "filename", "validation_summary", "status", "comment")
    list_filter = (('lab', RelatedFieldListFilter), ('status', AllValuesChoicesFieldListFilter))
    exclude = ('validation_list',)  # excludes validation_list that can be too big

    def is_readonly_field(self, f) -> bool:
        if f.name in ("url", "filename", "file_size"):
            return True
        if f.name == "user":
            return False
        return super().is_readonly_field(f)

    @admin_action("Process (Wait & Validate Only)")
    def process(self, request, queryset: QuerySet[UploadedClassificationsUnmapped]):
        for ufl in queryset:
            # noinspection PyArgumentList
            ClassificationImportMapInsertTask.run(upload_classifications_unmapped_id=ufl.pk, import_records=False)
            ufl.refresh_from_db()
            if validation_summary := ufl.validation_summary:
                if fatal_error := validation_summary.get("fatal_error"):
                    self.message_user(request, message=f"Fatal error: {fatal_error}", level=messages.ERROR)
                else:
                    records_mapped = validation_summary.get("records_mapped") or 0
                    self.message_user(request, message=f"{records_mapped} records mapped")

    @admin_action("Process (Async)")
    def process_async(self, request, queryset: QuerySet[UploadedClassificationsUnmapped]):
        for ufl in queryset:
            task = ClassificationImportMapInsertTask.si(ufl.pk)
            task.apply_async()


@admin.register(ResolvedVariantInfo)
class ResolvedVariantInfoAdmin(ModelAdminBasics):
    list_display = (
        'allele_info',
        'genome_build',
        'variant',
        'c_hgvs',
        'gene_symbol',
        'transcript_version',
        'genomic_sort',
        'error'
    )

    def has_add_permission(self, request):
        return False


class MatchingOnFilter(admin.SimpleListFilter):
    title = "Matching On"
    parameter_name = "field"

    def lookups(self, request, model_admin):
        return [("c_hgvs", "c.HGVS"), ("g_hgvs", "g.HGVS")]

    def queryset(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        if self.value() == "c_hgvs":
            queryset = queryset.filter(imported_c_hgvs__isnull=False)
        elif self.value() == "g_hgvs":
            queryset = queryset.filter(imported_g_hgvs__isnull=False)

        return queryset


class ImportedAlleleInfoValidationInline(admin.TabularInline):
    model = ImportedAlleleInfoValidation
    fields = ['c_hgvs_37', 'c_hgvs_38', 'confirmed', 'include', 'validation_tags']
    show_change_link = True

    def is_readonly_field(self, f):
        return True

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


class ImportedAlleleInfoDirtyFilter(admin.SimpleListFilter):
    title = "Dirt"
    parameter_name = "dirt"

    def lookups(self, request, model_admin):
        return [("dirty", "Dirty")]

    def queryset(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        if self.value() == "dirty":
            return queryset.filter(dirty_message__isnull=False)


class ImportedAlleleValidationFilter(admin.SimpleListFilter):
    title = 'Validations'
    parameter_name = 'validation_x'
    VALIDATION_TO_PATH = {
        "Missing 37": "builds__missing_37",
        "Missing 38": "builds__missing_38",
        "Transcript ID Change (normal)": "normalize__transcript_id_change",
        "Transcript ID Change (liftover)": "liftover__transcript_id_change",
        "Transcript Ver Change (normal)": "normalize__transcript_version_change",
        "Transcript Ver Change (liftover)": "liftover__transcript_version_change",
        "Gene Symbol Change (normal)": "normalize__gene_symbol_change",
        "Gene Symbol Change (liftover)": "liftover__gene_symbol_change",
        "General Issue (no resolve, unsup transcript)": "general"
    }

    default_value = None

    def lookups(self, request, model_admin):
        return list((item[1], item[0]) for item in ImportedAlleleValidationFilter.VALIDATION_TO_PATH.items())

    def queryset(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        if filter_path := self.value():
            return queryset.filter(**{f"latest_validation__validation_tags__{filter_path}__isnull": False})


@admin.register(ImportedAlleleInfo)
class ImportedAlleleInfoAdmin(ModelAdminBasics):
    list_display = (
        "id",
        "imported_hgvs",
        "imported_genome_build_patch_version",
        "status",
        "validation_include",
        "grch37_limit",
        "grch38_limit",
        "latest_validation",
        "variant_coordinate",
        "dirty_message"
        #
        # "created"
    )
    list_filter = (
        'imported_genome_build_patch_version', 'status', 'latest_validation__confirmed',
        ImportedAlleleInfoDirtyFilter,
        MatchingOnFilter,
        ImportedAlleleValidationFilter
    )
    search_fields = ('id', 'imported_c_hgvs', 'imported_g_hgvs', 'message')
    inlines = (ImportedAlleleInfoValidationInline,)

    @admin_list_column("grch37", "grch37")
    def grch37_limit(self, obj: ImportedAlleleInfo):
        return obj.grch37

    @admin_list_column("grch38", "grch38")
    def grch38_limit(self, obj: ImportedAlleleInfo):
        return obj.grch38

    @admin_list_column("Imported HGVS", "imported_c_hgvs")
    def imported_hgvs(self, obj: ImportedAlleleInfo):
        return obj.imported_c_hgvs or obj.imported_g_hgvs

    @admin_list_column("Include in Exports", "latest_validation__include", is_boolean=True)
    def validation_include(self, obj: ImportedAlleleInfo):
        if latest_validation := obj.latest_validation:
            return latest_validation.include
        return False

    # @admin_list_column("Latest Validation", "dirty_message")
    # def latest_validation_formatted(self, obj: ImportedAlleleInfo):
    #     if message := obj.latest_validation:
    #         escaped = escape(str(message))
    #         escaped.replace("\n", "<br/>")
    #         return SafeString(escaped)
    #
    # @admin_list_column("Dirty Message", "dirty_message")
    # def dirty_message_formatted(self, obj: ImportedAlleleInfo):
    #     if message := obj.dirty_message:
    #         escaped = escape(message)
    #         escaped.replace("\n", "<br/>")
    #         return SafeString(escaped)

    @admin_action("Dirty Check")
    def dirty_check(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        variant_matching_dry_run(queryset)

    @admin_action("Re-Match Soft")
    def re_match_soft(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        """
        Soft remwatch will leave everything linked while attempting to match again
        """
        for allele_info in queryset:
            allele_info.update_variant_coordinate()
            allele_info.refresh_and_save(force_update=True)
            allele_info.classification_import = None
            allele_info.status = ImportedAlleleInfoStatus.PROCESSING
            allele_info.save()

        rematched = reattempt_variant_matching(request.user, queryset, False)
        self.message_user(request, message=f"Allele Infos rematched {rematched}")

    @admin_action("Re-Match Hard (unmatched, rematches)")
    def re_match_hard(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        """
        Hard rematch will reset the imported allele info, and match from the start
        """
        rematched = reattempt_variant_matching(request.user, queryset, True)
        self.message_user(request, message=f"Allele Infos rematched {rematched}")

    @admin_action("Recalc c.hgvs")
    def update_variant_coordinate(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        for allele_info in queryset:
            allele_info.refresh_and_save(force_update=True)

    @admin_action("Validate")
    def validate(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        for allele_info in queryset:
            allele_info.apply_validation(force_update=True)
            allele_info.save()

    @admin_action("Confirm Include")
    def override_approve(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        queryset = queryset.select_related('latest_validation')
        iai: ImportedAlleleInfo
        for iai in queryset:
            iaiv = iai.latest_validation
            iaiv.include = True
            iaiv.confirmed = True
            iaiv.confirmed_by = request.user
            iaiv.save()

    @admin_action("Confirm Exclude")
    def override_reject(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        queryset = queryset.select_related('latest_validation')
        iai: ImportedAlleleInfo
        for iai in queryset:
            iaiv = iai.latest_validation
            iaiv.include = False
            iaiv.confirmed = True
            iaiv.confirmed_by = request.user
            iaiv.save()

    @admin_action("Confirm Reset to Default")
    def remove_confirmation(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        queryset = queryset.select_related('latest_validation')
        iai: ImportedAlleleInfo
        for iai in queryset:
            iavi = iai.latest_validation
            iavi.include = ImportedAlleleInfoValidation.should_include(iavi.validation_tags)
            iavi.confirmed = False
            iavi.confirmed_by = None
            iavi.save()

    @admin_action("Mark as Completed")
    def mark_as_completed(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        iai: ImportedAlleleInfo
        marked_complete = 0
        marked_failed = 0
        for iai in queryset.filter(
                status__in=(ImportedAlleleInfoStatus.PROCESSING, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD)):
            if iai.allele_id and (iai.grch37_id or iai.grch38_id):
                iai.status = ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS
                marked_complete += 1
            else:
                iai.status = ImportedAlleleInfoStatus.FAILED
                marked_failed += 1
            iai.save()
        self.message_user(request,
                          message=f'Records now marked as complete {marked_complete}, records marked as failure {marked_failed}',
                          level=messages.INFO)

    def has_add_permission(self, request):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(ImportedAlleleInfoValidation)
class ImportedAlleleInfoValidationAdmin(ModelAdminBasics):

    def is_readonly_field(self, f) -> bool:
        if f.name == 'confirmed_by_note':
            return False
        return True


class ClassificationGroupingEntryAdmin(admin.TabularInline):
    model = ClassificationGroupingEntry

    def has_add_permission(self, request, obj):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


class ClassificationGroupingSearchTermAdmin(admin.TabularInline):
    model = ClassificationGroupingSearchTerm

    def has_add_permission(self, request, obj):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


@admin.register(ClassificationGrouping)
class ClassificationGroupingAdmin(ModelAdminBasics):
    inlines = (ClassificationGroupingEntryAdmin, ClassificationGroupingSearchTermAdmin)
    list_display = ("pk", "classification_count", "allele", "lab", "allele_origin_bucket", "pathogenic_difference", "somatic_difference", "dirty")
    list_filter = ("lab", "allele_origin_bucket", "pathogenic_difference", "somatic_difference", "dirty")

    # @admin_list_column("gene_symbols")
    # def gene_symbols(self, obj: ClassificationGrouping):
    #     return ", ".join(obj.classificationgroupinggenesymbol_set.values_list("gene_symbol", flat=True))

    @admin_action("Refresh")
    def refresh(self, request, queryset: QuerySet[ClassificationGrouping]):
        queryset.update(dirty=True)
        for cg in queryset:
            cg.update()

    @admin_model_action(url_slug="refresh_all/", short_description="Refresh All", icon="fa-solid fa-arrows-rotate")
    def refresh_all(self, request):
        # ClassificationGrouping.objects.all().delete()
        for classification in Classification.objects.iterator():
            ClassificationGrouping.assign_grouping_for_classification(classification)

        ClassificationGrouping.objects.update(dirty=True)
        ClassificationGrouping.update_all_dirty()

    @admin_model_action(url_slug="refresh_all/", short_description="Refresh Dirty", icon="fa-solid fa-splotch")
    def refresh_dirty(self, request):
        ClassificationGrouping.update_all_dirty()


class ClassificationGroupingTabularAdmin(TabularInline):
    model = ClassificationGrouping

    def has_add_permission(self, request, obj):
        return False

    def has_delete_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(AlleleOriginGrouping)
class AlleleOriginGroupingAdmin(ModelAdminBasics):
    list_display = ("allele_grouping", "dirty")
    inlines = (ClassificationGroupingTabularAdmin,)

    @admin_model_action(url_slug="refresh_all/", short_description="Refresh All", icon="fa-solid fa-arrows-rotate")
    def refresh_all(self, request):
        AlleleOriginGrouping.objects.update(dirty=True)
        ClassificationGrouping.update_all_dirty()


class AlleleOriginGroupingTabularAdmin(TabularInline):
    model = AlleleOriginGrouping

    def has_add_permission(self, request, obj):
        return False

    def has_delete_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(AlleleGrouping)
class AlleleGroupingAdmin(ModelAdminBasics):
    inlines = (AlleleOriginGroupingTabularAdmin,)
    pass
