from datetime import timedelta

from django.contrib import admin
from django.contrib.admin import TabularInline
from django.db.models import QuerySet
from django.utils.safestring import SafeString

from annotation import models
from annotation.clinvar_xml_parser import ClinVarFetchRequest
from annotation.models import Citation, CitationFetchRequest, ClinVarRecordCollection, ClinVarRecord, VariantAnnotation, \
    ClinVar
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column, get_admin_url
from snpdb.models import VariantAllele

admin.site.register(models.AnnotationRun)
admin.site.register(models.AnnotationVersion)
admin.site.register(models.VariantAnnotationVersion)


@admin.register(ClinVar)
class ClinVarAdmin(ModelAdminBasics):

    list_display = ("pk", "version", "clinvar_variation_id", "variant", "clinvar_review_status")

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request):
        return False


class ClinVarRecordAdmin(TabularInline):
    model = ClinVarRecord

    fields = ("record_id", "submitter", "stars", "date_last_evaluated", "clinical_significance")

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request, obj=None):
        return False


@admin.register(ClinVarRecordCollection)
class ClinVarRecordCollectionAdmin(ModelAdminBasics):
    inlines = (ClinVarRecordAdmin, )
    list_per_page = 20

    # list_display = ("pk", "clinvar", "allele", "min_stars_loaded", "last_loaded")

    list_display = ("pk", "clinvar_variation_id", "min_stars_loaded", "last_loaded")

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

    @admin_list_column(limit=0)
    def allele(self, obj: ClinVarRecordCollection):
        try:
            allele = ClinVar.objects.filter(clinvar_variation_id=obj.clinvar_variation_id).order_by('-version').first().variant.allele
            href = get_admin_url(allele)
            return SafeString(f"<a href=\"{href}\">{allele}</a>")
        except Exception as ex:
            return str(ex)
    """

    def has_change_permission(self, request, obj=None):
        return False

    def has_add_permission(self, request):
        return False

    @admin_action("Refresh: If Old (current stars)")
    def refresh_old(self, request, queryset: QuerySet[ClinVarRecordCollection]):
        for obj in queryset:
            ClinVarFetchRequest(
                clinvar_variation_id=obj.clinvar_variation_id,
                min_stars=obj.min_stars_loaded
            ).fetch()

    @admin_action("Refresh: Force (current stars)")
    def refresh_force(self, request, queryset: QuerySet[ClinVarRecordCollection]):
        for obj in queryset:
            ClinVarFetchRequest(
                clinvar_variation_id=obj.clinvar_variation_id,
                min_stars=obj.min_stars_loaded,
                max_cache_age=timedelta(seconds=0)
            ).fetch()

    @admin_action("Refresh: Force (1+ stars)")
    def refresh_force_one_plus_stars(self, request, queryset: QuerySet[ClinVarRecordCollection]):
        for obj in queryset:
            ClinVarFetchRequest(
                clinvar_variation_id=obj.clinvar_variation_id,
                min_stars=1,
                max_cache_age=timedelta(seconds=0)
            ).fetch()


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
