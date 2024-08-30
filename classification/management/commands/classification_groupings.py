from django.core.management import BaseCommand

from classification.models import Classification
from classification.models.classification_grouping import ClassificationGrouping


class Command(BaseCommand):

    def handle(self, *args, **options):

        for classification in Classification.objects.iterator():
            ClassificationGrouping.assign_grouping_for_classification(classification)

        ClassificationGrouping.objects.all().update(dirty=True)
        ClassificationGrouping.update_all_dirty()