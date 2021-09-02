import json
from typing import Set

from django.contrib import admin, messages
from django.contrib.admin import RelatedFieldListFilter
from django.db.models import QuerySet
from django.utils import timezone

from annotation.models.models import AnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.classification_import import reattempt_variant_matching
from classification.enums.classification_enums import EvidenceCategory, SpecialEKeys, SubmissionSource, ShareLevel
from classification.models import EvidenceKey, EvidenceKeyMap, DiscordanceReport, DiscordanceReportClassification, \
    send_discordance_notification, ClinicalContext, ClassificationReportTemplate, ClassificationModification, \
    UploadedFileLab
from classification.models.classification import Classification
from library.guardian_utils import admin_bot
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column
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


@admin.register(Classification)
class ClassificationAdmin(ModelAdminBasics):
    list_display = ['id', 'lab', 'lab_record_id', 'share_level', 'clinical_significance', 'clinical_context', 'imported_genome_build', 'imported_c_hgvs', 'chgvs_grch37', 'chgvs_grch38', 'withdrawn', 'user', 'created_detailed', 'modified_detailed']
    list_filter = (('lab__organization', RelatedFieldListFilter), ('lab', RelatedFieldListFilter), ClassificationShareLevelFilter, VariantMatchedFilter, ClinicalContextFilter, ClassificationImportedGenomeBuildFilter, ('user', RelatedFieldListFilter),)
    search_fields = ('id', 'lab_record_id')
    list_per_page = 500
    inlines = (ClassificationModificationAdmin,)

    @admin_list_column(short_description="Created", order_field="created")
    def created_detailed(self, obj: Classification):
        return obj.created.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

    @admin_list_column(short_description="Modified", order_field="modified")
    def modified_detailed(self, obj: Classification):
        return obj.modified.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

    def has_add_permission(self, request):
        return False

    @admin_action("Populate base variant data")
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

    @admin_action("Revalidate")
    def revalidate(self, request, queryset):
        for vc in queryset:
            vc.revalidate(request.user)
        self.message_user(request, str(queryset.count()) + " records revalidated")

    @admin_action("Publish - Logged in Users")
    def publish_logged_in_users(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.ALL_USERS)

    @admin_action("Publish - Organisation")
    def publish_org(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.INSTITUTION)

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

    @admin_action("Variant re-matching")
    def reattempt_variant_matching(self, request, queryset: QuerySet[Classification]):
        valid_record_count, invalid_record_count = reattempt_variant_matching(request.user, queryset)
        if invalid_record_count:
            self.message_user(request, f'Records with missing or invalid builds/coordinates : {invalid_record_count}')
        self.message_user(request, f'Records revalidating : {valid_record_count}')

    @admin_action("Re-calculate cached chgvs")
    def recalculate_cached_chgvs(self, request, queryset: QuerySet[Classification]):
        for vc in queryset:
            vc.update_cached_c_hgvs()
            vc.save()

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

    @admin_action("Withdraw")
    def withdraw_true(self, request, queryset: QuerySet[Classification]):
        count = self.set_withdraw(request, queryset, True)
        self.message_user(request, f"{count} records now newly set to withdrawn")

    @admin_action("Un-Withdraw")
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
        return super(ClassificationAdmin, self).get_form(request, obj, widgets={
            'lab_record_id': admin.widgets.AdminTextInputWidget(),
            'chgvs_grch37': admin.widgets.AdminTextInputWidget(),
            'chgvs_grch37_full': admin.widgets.AdminTextInputWidget(),
            'chgvs_grch38': admin.widgets.AdminTextInputWidget(),
            'chgvs_grch38_full': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


@admin.register(ClinicalContext)
class ClinicalContextAdmin(ModelAdminBasics):
    list_display = ('id', 'allele', 'name', 'status', 'modified',)

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget()
        })

    def has_add_permission(self, request):
        return False

    @admin_action("Recalculate Status")
    def recalculate(self, request, queryset):
        for dc in queryset:
            dc.recalc_and_save(cause='Admin recalculation')  # cause of None should change to Unknown, which is accurate if this was required
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
        pretty = [e_key_cs.pretty_value(cs) for cs in sorted_list]
        return pretty

    @admin_list_column()
    def days_open(self, obj: DiscordanceReport) -> int:
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
        labs = sorted(labs, key=lambda x: x.name)
        return ", ".join([lab.name for lab in labs])

    def has_add_permission(self, request):
        return False

    @admin_action("Send discordance notification")
    def send_notification(self, request, queryset):
        ds: DiscordanceReport
        for ds in queryset:
            send_discordance_notification(ds)


@admin.register(UploadedFileLab)
class UploadedFileLabAdmin(ModelAdminBasics):
    list_display = ("pk", "lab", "created", "url", "status", "comment")
    pass

    def is_readonly_field(self, f) -> bool:
        if f.name in ("url", "filename"):
            return True
        return super().is_readonly_field(f)