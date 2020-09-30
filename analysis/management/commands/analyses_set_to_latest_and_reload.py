from django.core.management.base import BaseCommand
from django.db.models import F

from analysis.models import NodeStatus
from analysis.models.models_analysis import Analysis
from annotation.models.models import AnnotationVersion


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', help='Reload *everything*')

    def handle(self, *args, **options):
        qs = Analysis.objects.all()
        if not options["all"]:
            # During the VEP upgrade, VariantAnnotationVersion for IVAT was deleted, so these will be None
            qs = qs.filter(annotation_version__variant_annotation_version__isnull=True)

        for analysis in qs:
            print(f"Reloading: {analysis}..")
            analysis.annotation_version = AnnotationVersion.latest(analysis.genome_build)
            analysis.save()

            analysis.analysisnode_set.update(status=NodeStatus.DIRTY,
                                             count=None,
                                             version=F("version") + 1)

        print("Please run analyses_process_incomplete_nodes to actually do the reload")
