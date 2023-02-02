from django.core.management.base import BaseCommand

from analysis.models import AnalysisNode, NodeCache


class Command(BaseCommand):
    """
        @See https://github.com/SACGF/variantgrid/issues/696
        We can remove these now (were keeping them around in case we needed to roll back)
    """
    def handle(self, *args, **options):

        if total_nc := NodeCache.objects.count():
            print(f"Deleting legacy node cache, approx size = {total_nc}")
            for i, nc in enumerate(NodeCache.objects.all().select_related("node_version")):
                node = AnalysisNode.objects.get_subclass(pk=nc.node_version.node_id)
                # Should be false for everything else due to settings.ANALYSIS_NODE_CACHE_DB=False
                if not node.use_cache:
                    nc.delete()

                if i and i % 1000 == 0:
                    print(f"processed {i}/{total_nc} ({100 * i / total_nc}%)")
