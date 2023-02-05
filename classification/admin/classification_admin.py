import json
from typing import Set, Union

from django.contrib import admin, messages
from django.contrib.admin import RelatedFieldListFilter, BooleanFieldListFilter
from django.db.models import QuerySet, Q
from django.utils import timezone

from annotation.models.models import AnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.classification_import import reattempt_variant_matching
from classification.enums.classification_enums import EvidenceCategory, SpecialEKeys, SubmissionSource, ShareLevel
from classification.models import EvidenceKey, EvidenceKeyMap, DiscordanceReport, DiscordanceReportClassification, \
    send_discordance_notification, ClinicalContext, ClassificationReportTemplate, ClassificationModification, \
    UploadedClassificationsUnmapped, ClinicalContextRecalcTrigger, ImportedAlleleInfo, ClassificationImport
from classification.models.classification import Classification
from classification.models.classification_import_run import ClassificationImportRun, ClassificationImportRunStatus
from classification.models.classification_variant_info_models import ResolvedVariantInfo, ImportedAlleleInfoValidation
from classification.tasks.classification_import_map_and_insert_task import ClassificationImportMapInsertTask
from library.guardian_utils import admin_bot
from library.log_utils import log_admin_change
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column, AllValuesChoicesFieldListFilter, \
    admin_model_action
from snpdb.models import GenomeBuild, Lab


class VariantMatchedFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'Variant Match Status Filter'
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
    title = 'Share Level Filter'
    parameter_name = 'share_level'
    default_value = None

    def lookups(self, request, model_admin):
        return [(sl.key, sl.label) for sl in ShareLevel.ALL_LEVELS]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(share_level=self.value())
        return queryset


class ClinicalContextFilter(admin.SimpleListFilter):
    title = 'Clinical Context Filter'
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
    title = 'Imported Genome Build Filter'
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


class ClassificationModificationAdmin(admin.TabularInline):
    model = ClassificationModification

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(ClassificationImportRun)
class ClassificationImportRunAdmin(ModelAdminBasics):
    # change_list_template = 'classification/admin/change_list.html'
    list_display = ['id', 'identifier', 'row_count', 'status', 'from_file', 'created_detailed', 'modified_detailed']
    list_filter = ('status', )

    @admin_model_action(url_slug="create_dummy/", short_description="Create Dummy Import", icon="fa-solid fa-plus")
    def create_dummy(self, request):
        ClassificationImportRun(identifier="Dummy Import").save()
        self.message_user(request, "Created a dummy import")

    @admin_model_action(url_slug="mark_unfinished/", short_description="Mark OnGoing Imports as Unfinished", icon="fa-regular fa-trash-can")
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
        ('withdrawn', BooleanFieldListFilter),
        ClassificationShareLevelFilter,
        VariantMatchedFilter,
        ClinicalContextFilter,
        ClassificationImportedGenomeBuildFilter,
        ('user', RelatedFieldListFilter),
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
                data = get_evidence_fields_for_variant(genome_build, vc.variant, refseq_transcript_id, ensembl_transcript_id,
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

    @admin_action("Fixes: Revalidate")
    def revalidate(self, request, queryset):
        for vc in queryset:
            vc.revalidate(request.user)
        self.message_user(request, str(queryset.count()) + " records revalidated")

    @admin_action("Matching: Re-Match Variant")
    def reattempt_variant_matching(self, request, queryset: QuerySet[Classification]):
        for classification in queryset:
            _, created = classification.ensure_allele_info_with_created()
            if created:
                classification.save()

        allele_info_ids = queryset.values_list('allele_info', flat=True)
        allele_info_qs = ImportedAlleleInfo.objects.filter(id__in=allele_info_ids)

        valid_record_count, invalid_record_count = reattempt_variant_matching(request.user, allele_info_qs)
        if invalid_record_count:
            self.message_user(request, f'Records with missing or invalid builds/coordinates : {invalid_record_count}')
        self.message_user(request, f'Records re-matching : {valid_record_count}')

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
            if vc.share_level_enum >= share_level:
                already_published += 1
            elif vc.has_errors():
                in_error += 1
            else:
                try:
                    vc.publish_latest(user=request.user, share_level=share_level)
                    published += 1
                except ValueError as ve:
                    self.message_user(request, message='Error when publishing - ' + str(ve),
                                      level=messages.ERROR)
        if already_published:
            self.message_user(request, message=f"({already_published}) records had been previously published", level=messages.INFO)
        if in_error:
            self.message_user(request, message=f"({in_error}) records can't be published due to validation errors", level=messages.ERROR)
        self.message_user(request, message=f"({published}) records have been freshly published", level=messages.INFO)

    @admin_action("Publish: Organisation")
    def publish_org(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.INSTITUTION)

    @admin_action("Publish: Logged in Users")
    def publish_logged_in_users(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.ALL_USERS)

    @admin_action("State: Make Mutable")
    def make_mutable(self, request, queryset: QuerySet[Classification]):
        for vc in queryset:
            vc.patch_value(patch={}, user=request.user, source=SubmissionSource.VARIANT_GRID, remove_api_immutable=True, save=True)

    def set_withdraw(self, request, queryset: QuerySet[Classification], withdraw: bool) -> int:
        count = 0
        for vc in queryset:
            try:
                actioned = vc.set_withdrawn(user=request.user, withdraw=withdraw)
                if actioned:
                    count += 1
            except BaseException:
                pass
        return count

    @admin_action("State: Withdraw")
    def withdraw_true(self, request, queryset: QuerySet[Classification]):
        count = self.set_withdraw(request, queryset, True)
        self.message_user(request, f"{count} records now newly set to withdrawn")

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
    list_display = ('id', 'allele', 'name', 'status', 'modified', 'pending_cause', 'pending_status')
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
            dc.recalc_and_save(cause='Admin recalculation', cause_code=ClinicalContextRecalcTrigger.ADMIN)  # cause of None should change to Unknown, which is accurate if this was required
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
    title = 'Max Share Level Filter'
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
    search_fields = ('key',)

    fieldsets = (
        ('Basic', {'fields': ('key', 'label', 'sub_label')}),
        ('Position', {'fields': ('hide', 'evidence_category', 'order')}),
        ('Type', {'fields': ('mandatory', 'value_type', 'options', 'allow_custom_values', 'default_crit_evaluation')}),
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
                good_count = good_count + 1
            except ValueError as ve:
                self.message_user(request, str(ve), 'error')

        self.message_user(request, str(good_count) + " keys passed validation")

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'key': admin.widgets.AdminTextInputWidget(),
            'label': admin.widgets.AdminTextInputWidget(),
            'sub_label': admin.widgets.AdminTextInputWidget(),
            'see': admin.widgets.AdminURLFieldWidget(),
            'if_key': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


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
    readonly_fields = ["classification_original", "clinical_context_effective", "clinical_context_final", "withdrawn_final"]
    fields = ["classification_original", "clinical_context_effective", "clinical_context_final", "withdrawn_final"]

    def clinical_context_effective(self, drc: DiscordanceReportClassification):
        return drc.clinical_context_effective

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


@admin.register(DiscordanceReport)
class DiscordanceReportAdmin(ModelAdminBasics):
    list_display = ["pk", "allele", "report_started_date", "days_open", "classification_count", "clinical_sigs", "labs"]
    list_select_related = ('clinical_context', 'clinical_context__allele')
    list_filter = [DiscordanceReportAdminLabFilter]
    inlines = (DiscordanceReportClassificationAdmin,)

    @admin_list_column("Allele", order_field="clinical_context__allele__pk")
    def allele(self, obj: DiscordanceReport) -> str:
        cc = obj.clinical_context
        return str(cc.allele)

    @admin_list_column("Clinical Significances")
    def clinical_sigs(self, obj: DiscordanceReport) -> str:
        clinical_sigs = set()
        for dr in DiscordanceReportClassification.objects.filter(report=obj):
            clinical_sigs.add(dr.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

        e_key_cs = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        sorter = e_key_cs.classification_sorter_value
        sorted_list = sorted(clinical_sigs, key=lambda x: (sorter(x), x))
        pretty = ", ".join([e_key_cs.pretty_value(cs) for cs in sorted_list])
        return pretty

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
        labs: Set[Lab] = set()
        drc: DiscordanceReportClassification
        for drc in obj.discordancereportclassification_set.all():
            if c := drc.classification_original.classification:
                labs.add(c.lab)
        labs_sorted = sorted(labs, key=lambda x: x.name)
        return ", ".join([lab.name for lab in labs_sorted])

    def has_add_permission(self, request):
        return False

    @admin_action("Send discordance notification")
    def send_notification(self, request, queryset):
        ds: DiscordanceReport
        for ds in queryset:
            send_discordance_notification(ds, cause="Admin Send discordance notification")

    @admin_action("Re-calculate latest")
    def re_calculate(self, request, queryset):
        ds: DiscordanceReport
        for ds in queryset:
            ds.clinical_context.recalc_and_save(cause="Admin recalculation", cause_code=ClinicalContextRecalcTrigger.ADMIN)


@admin.register(UploadedClassificationsUnmapped)
class UploadedClassificationsUnmappedAdmin(ModelAdminBasics):
    list_display = ("pk", "lab", "created", "filename", "validation_summary", "status", "comment")
    list_filter = (('lab', RelatedFieldListFilter), ('status', AllValuesChoicesFieldListFilter))

    def is_readonly_field(self, f) -> bool:
        if f.name in ("url", "filename"):
            return True
        return super().is_readonly_field(f)

    @admin_action("Process (Wait & Validate Only)")
    def process(self, request, queryset: QuerySet[UploadedClassificationsUnmapped]):
        for ufl in queryset:
            ClassificationImportMapInsertTask.run(upload_classifications_unmapped_id=ufl.pk, import_records=False)

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


class ValidationFilter(admin.SimpleListFilter):
    title = "Validation"
    parameter_name = "validation"

    def lookups(self, request, model_admin):
        return [("gene_symbol", "Gene Symbol")]

    def queryset(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        if self.value() == "gene_symbol":
            return queryset.filter(
                Q(latest_validation__validation_tags__normalize__gene_symbol_change__isnull=False) | \
                Q(latest_validation__validation_tags__liftover__gene_symbol_change__isnull=False)
            ).filter(
                Q(latest_validation__validation_tags__normalize__c_nomen_change__isnull=True) & \
                Q(latest_validation__validation_tags__liftover__c_nomen_change__isnull=True) & \
                Q(latest_validation__validation_tags__builds__missing_37__isnull=True) & \
                Q(latest_validation__validation_tags__builds__missing_38__isnull=True)
            )


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


@admin.register(ImportedAlleleInfo)
class ImportedAlleleInfoAdmin(ModelAdminBasics):
    list_display = (
        "id",
        "imported_hgvs",
        "imported_genome_build_patch_version",
        "status",
        "validation_include",
        "grch37",
        "grch38",
        "latest_validation"
        #"variant_coordinate",
        #"created"
    )
    list_filter = ('imported_genome_build_patch_version', 'status', 'latest_validation__confirmed', ValidationFilter, MatchingOnFilter)
    search_fields = ('imported_c_hgvs', 'imported_g_hgvs')
    inlines = (ImportedAlleleInfoValidationInline,)


    @admin_list_column("Imported HGVS", "imported_c_hgvs")
    def imported_hgvs(self, obj: ImportedAlleleInfo):
        return obj.imported_c_hgvs or obj.imported_g_hgvs

    @admin_list_column("Include in Exports", "latest_validation__include", is_boolean=True)
    def validation_include(self, obj: ImportedAlleleInfo):
        if latest_validation := obj.latest_validation:
            return latest_validation.include
        return False

    @admin_action("Re-Match")
    def re_match(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        reattempt_variant_matching(request.user, queryset)

        # classifications_qs = Classification.objects.filter(allele_info__in=queryset)
        # self.message_user(request, f'Found : {classifications_qs.count()} linked classifications')
        #valid_record_count, invalid_record_count = reattempt_variant_matching(request.user, classifications_qs)
        #if invalid_record_count:
        #    self.message_user(request, f'Records with missing or invalid builds/coordinates : {invalid_record_count}')
        #self.message_user(request, f'Classification Records re-matching : {valid_record_count}')

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

    @admin_action("Confirm Reset")
    def remove_confirmation(self, request, queryset: QuerySet[ImportedAlleleInfo]):
        queryset = queryset.select_related('latest_validation')
        iai: ImportedAlleleInfo
        for iai in queryset:
            iavi = iai.latest_validation
            iavi.include = ImportedAlleleInfoValidation.should_include(iavi.validation_tags)
            iavi.confirmed = False
            iavi.confirmed_by = None
            iavi.save()

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
