from time import sleep

from django.core.management.base import BaseCommand

from analysis.models import NodeStatus, AnalysisNode
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_utils import update_nodes


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--max-queue', type=int, default=200, help='Max number of queued nodes')

    def handle(self, *args, **options):
        max_queue = options["max_queue"]

        print("Assuming celery queues are EMPTY!")
        print("Marking loading nodes as dirty....")
        loading_nodes = AnalysisNode.objects.filter(status__in=NodeStatus.LOADING_STATUSES)
        loading_nodes.update(status=NodeStatus.DIRTY)
        dirty_nodes = AnalysisNode.objects.filter(status=NodeStatus.DIRTY)

        qs = Analysis.objects.filter(pk__in=dirty_nodes.values_list("analysis_id"))
        print(f"Reloading {dirty_nodes.count()} dirty nodes in analyses: {qs.count()} (max queue: {max_queue})")
        QUEUED_STATUS = [NodeStatus.QUEUED, NodeStatus.LOADING, NodeStatus.LOADING_CACHE]

        for analysis in qs:
            while True:
                num_queued = AnalysisNode.objects.filter(status__in=QUEUED_STATUS).count()
                if num_queued > max_queue:
                    print(f"waiting for queue {num_queued} to shrink below {max_queue}")
                    sleep(5)
                else:
                    break

            print(f"Updating nodes for {analysis}")
            update_nodes(analysis.pk)
