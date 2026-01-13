from django.core.management import BaseCommand

from classification.models import ClassificationGrouping, Overlap, OverlapContribution
from classification.services.overlaps_services import OverlapServices


class Command(BaseCommand):

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        Overlap.objects.all().delete()
        OverlapContribution.objects.all().delete()

        for cg in ClassificationGrouping.objects.iterator():
            OverlapServices.update_classification_grouping_overlap_contribution(cg)

        print(f"Overlap Count = {Overlap.objects.count()}")
        print(f"Overlap Contribution Count = {OverlapContribution.objects.count()}")
