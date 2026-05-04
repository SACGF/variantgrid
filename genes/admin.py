from django.contrib import admin, messages
from django.db.models import QuerySet, TextField
from django.db.models.functions import Cast
from django.utils.safestring import SafeString
from guardian.admin import GuardedModelAdmin

from genes import models
from genes.models import GeneSymbol, GeneCoverageCollection
from snpdb.admin_utils import ModelAdminBasics, admin_list_column, admin_action
from snpdb.archive import ArchivePreconditionError


@admin.register(models.GeneList)
class GeneListAdmin(GuardedModelAdmin):
    search_fields = ('name', )


@admin.register(models.GeneSymbol)
class GeneSymbolAdmin(ModelAdminBasics):
    search_fields = ('symbol_text', )
    list_display = ['symbol_text', 'genes']

    def get_queryset(self, request):
        qs: QuerySet[GeneSymbol] = super().get_queryset(request)
        qs = qs.annotate(symbol_text=Cast('symbol', output_field=TextField()))
        qs = qs.order_by('symbol_text')
        return qs

    @admin_list_column("symbol", order_field="symbol_text")
    def symbol_text(self, obj: GeneSymbol):
        return obj.symbol or SafeString("<i>BLANK</i>")

    @admin_list_column("genes")
    def genes(self, obj: GeneSymbol):
        return ", ".join(str(g) for g in obj.genes)

    @admin_action("panelapp au refresh")
    def panelapp_au_refresh(self, request, queryset: QuerySet[GeneSymbol]):
        from ontology.panel_app_ontology import update_gene_relations
        for gene_symbol in queryset:
            update_gene_relations(gene_symbol)


@admin.register(models.ActiveSampleGeneList)
class ActiveSampleGeneListAdmin(ModelAdminBasics):
    pass


@admin.register(models.GeneAnnotationRelease)
class GeneAnnotationReleaseAdmin(ModelAdminBasics):
    pass


@admin.register(models.GeneInfo)
class GeneInfoAdmin(ModelAdminBasics):
    pass


@admin.register(models.GeneListGeneSymbol)
class GeneListGeneSymbolAdmin(ModelAdminBasics):
    pass


@admin.register(models.GeneSymbolAlias)
class GeneSymbolAliasAdmin(ModelAdminBasics):
    pass


@admin.register(models.PanelAppPanel)
class PanelAppPanelAdmin(ModelAdminBasics):
    pass


@admin.register(models.PanelAppPanelLocalCache)
class PanelAppPanelLocalCacheAdmin(ModelAdminBasics):
    pass


@admin.register(models.SampleGeneList)
class SampleGeneListAdmin(ModelAdminBasics):
    pass


@admin.register(models.GeneCoverageCollection)
class GeneCoverageCollectionAdmin(ModelAdminBasics):
    list_display = ("pk", "path", "data_state", "genome_build", "data_archived_date")
    list_filter = ("data_state", "genome_build", "data_archived_date")
    actions = ["archive_selected", "restore_selected"]

    @admin_action("Archive selected")
    def archive_selected(self, request, queryset: QuerySet[GeneCoverageCollection]):
        from genes.archive import archive_gene_coverage_collection
        for gcc in queryset:
            try:
                archive_gene_coverage_collection(gcc, request.user, reason="Admin action")
                messages.info(request, f"Archived GCC pk={gcc.pk}")
            except ArchivePreconditionError as e:
                messages.error(request, f"GCC pk={gcc.pk}: {e}")

    @admin_action("Restore selected")
    def restore_selected(self, request, queryset: QuerySet[GeneCoverageCollection]):
        from genes.archive import restore_gene_coverage_collection
        for gcc in queryset:
            try:
                restore_gene_coverage_collection(gcc, request.user)
                messages.info(request, f"Restore queued for GCC pk={gcc.pk}")
            except ValueError as e:
                messages.error(request, f"GCC pk={gcc.pk}: {e}")
