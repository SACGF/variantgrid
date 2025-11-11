import argparse

from django.core.management import BaseCommand

from classification.models import Classification, ClassificationModification, ClassificationSummaryCalculator
from classification.models.classification_grouping import ClassificationGrouping, AlleleOriginGrouping


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--summary', required=False, action="store_true", help="Refreshes the summary data assigned to each classification")
        parser.add_argument('--refresh', required=False, action="store_true", help="Refreshes all existing groups, but not which classifications belong to them")
        parser.add_argument('--all', required=False, action="store_true", help="Refreshes which classification belongs to which group, and the groups, may take a long time")
        parser.add_argument('--dirty', required=False, action="store_true", help="Updates all records left in a dirty state")

    def handle(self, *args, **options):
        summary = options.get("summary")
        all = options.get("all")
        dirty = options.get("dirty")
        refresh = options.get("refresh")

        if not any((summary, all, dirty, refresh)):
            raise ValueError("Must provide one or more of summary, all, dirty, refresh")

        if summary:
            for index, cm in enumerate(ClassificationModification.objects.filter(is_last_published=True).select_related("classification").iterator()):
                classification = cm.classification
                classification.summary = ClassificationSummaryCalculator(cm).cache_dict()
                classification.save(update_fields=["summary"])
                if index % 1000 == 0 and index:
                    print(f"Updating {index} classification summaries")

        if all:
            for index, classification in enumerate(Classification.objects.iterator()):
                ClassificationGrouping.assign_grouping_for_classification(classification)
                if index % 1000 == 0 and index:
                    print(f"Updating {index} classification assigned to classification groups")
            print("About to update all classification groups")
            ClassificationGrouping.objects.all().update(dirty=True)

        if all or dirty or refresh:
            qs = ClassificationGrouping.objects.all()
            if not refresh:
                qs = qs.filter(dirty=True)

            index = 0
            for index, dirty in enumerate(qs.iterator()):
                dirty.update()
                if index % 1000 == 0 and index:
                    print(f"Updating {index} classification groupings")
            print(f"Updated {index+1} classification groupings")

            for index, dirty in enumerate(AlleleOriginGrouping.objects.filter(dirty=True).iterator()):
                dirty.update()
                if index % 1000 == 0 and index:
                    print(f"Updating {index} allele groupings")
            print(f"Updated {index+1} allele groupings")
