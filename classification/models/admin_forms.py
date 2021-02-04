from django.contrib import admin, messages
from django.contrib.auth.models import User
from django_admin_json_editor.admin import JSONEditorWidget
import json
from django.conf import settings

from annotation.models.models import AnnotationVersion
from library.guardian_utils import admin_bot
from snpdb.admin import ModelAdminBasics
from snpdb.models import ImportSource, Lab, Organization, GenomeBuild
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant
from classification.enums.classification_enums import EvidenceCategory, SpecialEKeys, SubmissionSource, ShareLevel
from classification.models import PatchMeta, EvidenceKey, email_discordance_for_classification, ConditionText, \
    ConditionTextMatch
from classification.models.classification import Classification, ClassificationImport
from classification.models.classification_patcher import patch_merge_age_units, patch_fuzzy_age
from classification.classification_import import process_classification_import


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


class ClassificationLabFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'Lab Filter'
    parameter_name = 'lab'
    default_value = None

    def lookups(self, request, model_admin):
        return [(lab.id, lab.name) for lab in Lab.objects.all()]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(lab=self.value())
        return queryset


class ClassificationOrgFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'Org Filter'
    parameter_name = 'org'
    default_value = None

    def lookups(self, request, model_admin):
        return [(org.id, org.name) for org in Organization.objects.all()]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(lab__organization=self.value())
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


class ClassificationUserFilter(admin.SimpleListFilter):
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


class ClassificationAdmin(admin.ModelAdmin):
    list_display = ['id', 'lab', 'lab_record_id', 'share_level', 'clinical_significance', 'clinical_context', 'imported_genome_build', 'withdrawn', 'chgvs_grch37', 'chgvs_grch38', 'user', 'created']
    list_filter = (ClassificationOrgFilter, ClassificationLabFilter, ClassificationShareLevelFilter, VariantMatchedFilter, ClinicalContextFilter, ClassificationImportedGenomeBuildFilter, ClassificationUserFilter,)
    search_fields = ('id', 'lab_record_id')

    fieldsets = (
        ('ID', {'fields': ('lab_record_id',)}),
        ('Links', {'fields': ('sample', 'lab', 'user', 'clinical_context')}),
        ('Variant', {'fields': ('chgvs_grch37', 'chgvs_grch38')}),
        ('Data', {'fields': ('clinical_significance', 'evidence',)}),
    )

    def populate_base_variant_data(self, request, queryset):
        for vc in queryset:
            refseq_transcript_id = vc.get(SpecialEKeys.REFSEQ_TRANSCRIPT_ID)
            ensembl_transcript_id = vc.get(SpecialEKeys.ENSEMBL_TRANSCRIPT_ID)
            if not vc.variant:
                self.message_user(request, "(%i) No variant for classification" % vc.id)
            else:
                genome_build = vc.variant.genome_build
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

    populate_base_variant_data.short_description = "Populate base variant data"

    def migration_fixes(self, request, queryset):
        for vc in queryset:
            vc.fix_migration_stuff(request.user)

        self.message_user(request, f"{queryset.count()} records migrated")

    migration_fixes.short_description = "Migration fixes"

    def age_fixes(self, request, queryset):

        def age_patches(patch: PatchMeta):
            patch_merge_age_units(patch)
            if settings.VARIANT_CLASSIFICATION_AUTOFUZZ_AGE:
                patch_fuzzy_age(patch)

        vc: Classification
        for vc in queryset:
            vc.patch_history(age_patches)

        self.message_user(request, f"{queryset.count()} records age changes")

    age_fixes.short_description = "** Age fixes"

    def revalidate(self, request, queryset):
        for vc in queryset:
            vc.revalidate(request.user)
        self.message_user(request, str(queryset.count()) + " records revalidated")

    revalidate.short_description = "Revalidate"

    def publish_logged_in_users(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.ALL_USERS)

    publish_logged_in_users.short_description = 'Publish - App Users'

    def publish_3rd_party_dbs(self, request, queryset):
        self.publish_share_level(request, queryset, ShareLevel.PUBLIC)

    publish_3rd_party_dbs.short_description = 'Publish - 3rd Party Databases'

    def publish_share_level(self, request, queryset, share_level: ShareLevel):
        already_published = 0
        in_error = 0
        published = 0
        vc: Classification
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

    def reattempt_variant_matching(self, request, queryset):
        qs = queryset
        qs.order_by('evidence__genome_build')

        invalid_genome_build_count = 0
        valid_record_count = 0
        imports_by_genome = {}

        for vc in qs:
            try:
                genome_build = vc.get_genome_build()
                if not genome_build.pk in imports_by_genome:
                    imports_by_genome[genome_build.pk] = ClassificationImport.objects.create(user=request.user,
                                                                                             genome_build=genome_build)
                vc_import = imports_by_genome[genome_build.pk]
                vc.set_variant(variant=None, message='Admin has re-triggered variant matching')
                vc.classification_import = vc_import
                vc.save()
                valid_record_count = valid_record_count + 1

            except:
                invalid_genome_build_count = invalid_genome_build_count + 1

        for vc_import in imports_by_genome.values():
            process_classification_import(vc_import, ImportSource.API)

        if invalid_genome_build_count:
            self.message_user(request, f'Records with missing or invalid genome_builds : {invalid_genome_build_count}')
        self.message_user(request, f'Records revalidating : {valid_record_count}')

    reattempt_variant_matching.short_description = 'Variant re-matching'

    def recalculate_cached_chgvs(self, request, queryset):
        vc: Classification
        for vc in queryset:
            vc.update_cached_c_hgvs()
            vc.save()

    recalculate_cached_chgvs.short_description = 'Re-calculate cached chgvs'

    def email_discordance_notification(self, request, queryset):
        emails_sent = False
        for vc in queryset:
            this_vc_sent = email_discordance_for_classification(vc)
            if this_vc_sent:
                emails_sent = True
        if emails_sent:
            self.message_user(request, 'Some emails were sent')
        else:
            self.message_user(request, 'No emails were sent')

    email_discordance_notification.short_description = 'Send Discordance Notifications'

    def set_withdraw(self, request, queryset, withdraw: bool) -> int:
        vc: Classification
        count = 0
        for vc in queryset:
            try:
                actioned = vc.set_withdrawn(user=request.user, withdraw=withdraw)
                if actioned:
                    count += 1
            except:
                pass
        return count

    def withdraw_true(self, request, queryset):
        count = self.set_withdraw(request, queryset, True)
        self.message_user(request, f"{count} records now newly set to withdrawn")

    withdraw_true.short_description = 'Withdraw'

    def withdraw_false(self, request, queryset):
        count = self.set_withdraw(request, queryset, False)
        self.message_user(request, f"{count} records now newly set to un-withdrawn")

    withdraw_false.short_description = 'Un-Withdraw'

    actions = [migration_fixes,
               revalidate,
               populate_base_variant_data,
               publish_logged_in_users,
               publish_3rd_party_dbs,
               reattempt_variant_matching,
               recalculate_cached_chgvs,
               email_discordance_notification,
               age_fixes,
               withdraw_true,
               withdraw_false]

    def get_form(self, request, obj=None, **kwargs):
        return super(ClassificationAdmin, self).get_form(request, obj, widgets={
            'lab_record_id': admin.widgets.AdminTextInputWidget()
        }, **kwargs)


class ClinicalContextAdmin(admin.ModelAdmin):
    list_display = ('id', 'allele', 'name', 'status', 'modified',)

    def recalc(self, request, queryset):
        for dc in queryset:
            dc.recalc_and_save(cause='Admin recalculation')  # cause of None should change to Unknown, which is accurate if this was required
        self.message_user(request, 'Recalced %i statuses' % queryset.count())

    recalc.short_description = 'Recalc status'
    actions = [recalc]


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


class EvidenceKeyAdmin(admin.ModelAdmin):
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

    validate.short_description = "Validate key"

    actions = [validate]

    def get_form(self, request, obj=None, **kwargs):

        options_schema = {
            'type': 'array',
            'title': 'select',
            'description': 'Options - ONLY meaningful if value type is Select',
            'items': {
                'type': 'object',
                'title': 'option',
                'properties': {
                    'key': {
                        'type': 'string',
                        'description': 'value stored against the record'
                    },
                    'label': {
                        'type': 'string',
                        'description': 'value displayed to the user'
                    },
                    'index': {
                        'type': 'number',
                        'required': True,
                        'description': 'Number to be used for RedCap exports'
                    }
                }
            }
        }
        examples_schema = {
            'type': 'array',
            'title': 'examples',
            'items': {
                'type': 'string',
                'title': 'example'
            }
        }

        return super(EvidenceKeyAdmin, self).get_form(request, obj, widgets={

            'key': admin.widgets.AdminTextInputWidget(),
            'label': admin.widgets.AdminTextInputWidget(),
            'sub_label': admin.widgets.AdminTextInputWidget(),
            'see': admin.widgets.AdminURLFieldWidget(),
            'if_key': admin.widgets.AdminTextInputWidget(),
            'options': JSONEditorWidget(options_schema, True),
            'examples': JSONEditorWidget(examples_schema, True),

        }, **kwargs)


class ClassificationReportTemplateAdmin(admin.ModelAdmin):
    list_display = ('name', 'modified')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'template': admin.widgets.AdminTextareaWidget()
        }, **kwargs)


class ConditionTextStatusFilter(admin.SimpleListFilter):
    title = 'Status Filter'
    parameter_name = 'status_filter'
    default_value = None

    def lookups(self, request, model_admin):
        return [("outstanding", "Only with outstanding classifications")]

    def queryset(self, request, queryset):
        if value := self.value():
            if value == "outstanding":
                queryset = queryset.exclude(classifications_count_outstanding=0)
        return queryset


class ConditionTextAdmin(ModelAdminBasics):
    search_fields = ('id', 'normalized_text')
    list_display = ["pk", "lab", "normalized_text", "classifications_count", "classifications_count_outstanding", "min_auto_match_score"]
    list_filter = [ConditionTextStatusFilter, ClassificationLabFilter]

    def auto_match(self, request, queryset):
        condition_text: ConditionText
        for condition_text in queryset:
            ConditionTextMatch.attempt_automatch(condition_text=condition_text, force=True, server_search=True)
            self.message_user(request, message=f"Automatching of {condition_text.normalized_text} resulted in min score of {condition_text.min_auto_match_score} matching {condition_text.classification_match_count}",
                              level=messages.INFO)

    auto_match.short_description = "Automatch (overwrite existing data)"
    actions = [auto_match]


class ConditionTextMatchAdmin(ModelAdminBasics):
    list_display = ["pk", "condition_text", "gene_symbol", "classification", "condition_xrefs", "condition_multi_operation", "created"]


class ClinVarExportAdmin(ModelAdminBasics):
    list_display = ["pk", "lab", "allele", "transcript", "gene_symbol", "created"]



class ConditionTextMatchAdmin(ModelAdminBasics):
    list_display = ["pk", "condition_text", "gene_symbol", "classification", "condition_xrefs", "condition_multi_operation", "created"]


class ClinVarExportAdmin(ModelAdminBasics):
    list_display = ["pk", "lab", "allele", "transcript", "gene_symbol", "created"]
