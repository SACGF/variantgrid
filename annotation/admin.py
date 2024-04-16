from datetime import timedelta, datetime

from django.contrib import admin, messages
from django.contrib.admin import TabularInline
from django.db.models import QuerySet
from django.utils.safestring import SafeString

from annotation import models
from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.clinvar_xml_parser import CLINVAR_RECORD_CACHE_DAYS
from annotation.clinvar_xml_parser_via_vcv import ClinVarXmlParserViaVCV
from annotation.models import Citation, CitationFetchRequest, ClinVarRecordCollection, ClinVarRecord, ClinVar, \
    AnnotationRun
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column, get_admin_url

admin.site.register(models.AnnotationVersion)
admin.site.register(models.GeneAnnotationVersion)
admin.site.register(models.VariantAnnotationVersion)


@admin.register(AnnotationRun)
class AnnotationRunAdmin(ModelAdminBasics):
    list_display = ("pk", "pipeline_type", "status")
    list_filter = ("pipeline_type", "status")


@admin.register(ClinVar)
class ClinVarAdmin(ModelAdminBasics):

    list_display = ("pk", "version", "clinvar_variation_id", "variant", "clinvar_review_status")

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request):
        return False


class ClinVarRecordAdmin(TabularInline):
    model = ClinVarRecord

    fields = ("record_id", "submitter", "hgvs", "condition", "stars", "date_last_evaluated", "clinical_significance", "somatic_clinical_significance", "allele_origin_bucket")
    readonly_fields = ("hgvs", )

    def hgvs(self, obj: ClinVarRecord):
        if obj.c_hgvs:
            return obj.c_hgvs
        return obj.variant_coordinate

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request, obj=None):
        return False


@admin.register(ClinVarRecordCollection)
class ClinVarRecordCollectionAdmin(ModelAdminBasics):
    search_fields = ("clinvar_variation_id", "allele__id")
    inlines = (ClinVarRecordAdmin, )
    list_per_page = 20

    list_display = ("pk", "clinvar_variation_id", "allele_link", "max_stars", "last_loaded")

    """
    # these took prohibitively long to load
    
    @admin_list_column(limit=0)
    def clinvar(self, obj: ClinVarRecordCollection):
        try:
            clinvar = ClinVar.objects.filter(clinvar_variation_id=obj.clinvar_variation_id).order_by('-version').first()
            href = get_admin_url(clinvar)
            return SafeString(f"<a href=\"{href}\">{clinvar.clinvar_variation_id}</a>")
        except Exception as ex:
            return str(ex)
    """

    @admin_list_column("allele", limit=0)
    def allele_link(self, obj: ClinVarRecordCollection):
        if allele := obj.allele:
            href = get_admin_url(allele)
            return SafeString(f"<a href=\"{href}\">{allele}</a>")

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request):
        return False

    @admin_action(f"Refresh: If Older than {CLINVAR_RECORD_CACHE_DAYS} days")
    def refresh_old(self, request, queryset: QuerySet[ClinVarRecordCollection]):
        for obj in queryset:
            ClinVarFetchRequest(
                clinvar_variation_id=obj.clinvar_variation_id,
            ).fetch()

    @admin_action("Refresh: Force (using VCV - default)")
    def refresh_force_vcv(self, request, queryset: QuerySet[ClinVarRecordCollection]):
        start = datetime.now()
        for obj in queryset:
            ClinVarFetchRequest(
                clinvar_variation_id=obj.clinvar_variation_id,
                max_cache_age=timedelta(seconds=0),
                parser=ClinVarXmlParserViaVCV
            ).fetch()
        duration = datetime.now() - start
        messages.info(request, message=f"Fetching took {duration}")

    # @admin_action("Refresh: Force (using RCVs)")
    # def refresh_force_rcvs(self, request, queryset: QuerySet[ClinVarRecordCollection]):
    #     start = datetime.now()
    #     for obj in queryset:
    #         ClinVarFetchRequest(
    #             clinvar_variation_id=obj.clinvar_variation_id,
    #             max_cache_age=timedelta(seconds=0),
    #             parser=ClinVarXmlParserViaRCVs
    #         ).fetch()
    #     duration = datetime.now() - start
    #     messages.info(request, message=f"Fetching took {duration}")


class HasErrorFilter(admin.SimpleListFilter):
    title = "Has Errors"
    parameter_name = "error_mode"

    def lookups(self, request, model_admin):
        return [("has_error", "Has Errors")]

    def queryset(self, request, queryset: QuerySet[Citation]):
        if self.value() == "has_error":
            queryset = queryset.filter(error__isnull=False)
        return queryset


class HasLoaded(admin.SimpleListFilter):
    title = "Has Loaded"
    parameter_name = "loaded"

    def lookups(self, request, model_admin):
        return [
            ("not-loaded", "Not-Loaded"),
            ("loaded", "Loaded")
        ]

    def queryset(self, request, queryset: QuerySet[Citation]):
        if self.value() == "not-loaded":
            queryset = queryset.filter(last_loaded__isnull=True)
        if self.value() == "loaded":
            queryset = queryset.filter(last_loaded__isnull=False)
        return queryset


@admin.register(Citation)
class CitationAdmin(ModelAdminBasics):
    list_display = ('id', 'title', 'error', 'last_loaded', 'old_id')
    list_filter = ('source', HasErrorFilter, HasLoaded)
    search_fields = ('id', 'source')

    @admin_action("Force Refresh")
    def force_refresh(self, request, queryset: QuerySet[Citation]):
        CitationFetchRequest.fetch_all_now(queryset, cache_age=timedelta(seconds=0))

    @admin_action("Load (if Unloaded)")
    def load_if_unloaded(self, request, queryset: QuerySet[Citation]):
        CitationFetchRequest.fetch_all_now(queryset)
