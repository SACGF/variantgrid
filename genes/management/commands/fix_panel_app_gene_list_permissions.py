from django.core.management import BaseCommand

from genes.models import PanelAppPanelLocalCacheGeneList
from library.guardian_utils import add_public_group_read_permission


class Command(BaseCommand):
    def handle(self, *args, **options):
        for panel_app_panel_local_cache_gene_list in PanelAppPanelLocalCacheGeneList.objects.all():
            gene_list = panel_app_panel_local_cache_gene_list.gene_list
            add_public_group_read_permission(gene_list)
