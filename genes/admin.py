from django.contrib import admin
from guardian.admin import GuardedModelAdmin

from genes import models


@admin.register(models.GeneList)
class GeneListAdmin(GuardedModelAdmin):
    search_fields = ('name', )


admin.site.register(models.ActiveSampleGeneList)
admin.site.register(models.GeneAnnotationRelease)
admin.site.register(models.GeneInfo)
admin.site.register(models.GeneListGeneSymbol)
admin.site.register(models.GeneSymbol)
admin.site.register(models.GeneSymbolAlias)
admin.site.register(models.PanelAppPanel)
admin.site.register(models.PanelAppPanelLocalCache)
admin.site.register(models.SampleGeneList)
