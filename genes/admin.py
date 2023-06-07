from django.contrib import admin
from django.db.models import QuerySet, TextField
from django.db.models.functions import Cast
from django.utils.safestring import SafeString
from guardian.admin import GuardedModelAdmin

from genes import models
from genes.models import GeneSymbol
from snpdb.admin_utils import ModelAdminBasics, admin_list_column


@admin.register(models.GeneList)
class GeneListAdmin(GuardedModelAdmin):
    search_fields = ('name', )


@admin.register(models.GeneSymbol)
class GeneSymbolAdmin(ModelAdminBasics):
    search_fields = ('name', )
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



admin.site.register(models.ActiveSampleGeneList)
admin.site.register(models.GeneAnnotationRelease)
admin.site.register(models.GeneInfo)
admin.site.register(models.GeneListGeneSymbol)
admin.site.register(models.GeneSymbolAlias)
admin.site.register(models.PanelAppPanel)
admin.site.register(models.PanelAppPanelLocalCache)
admin.site.register(models.SampleGeneList)
