from django.contrib import admin
from django.db.models import QuerySet, TextField
from django.db.models.functions import Cast
from django.utils.safestring import SafeString
from guardian.admin import GuardedModelAdmin

from genes import models
from genes.models import GeneSymbol
from snpdb.admin_utils import ModelAdminBasics, admin_list_column, admin_action


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
