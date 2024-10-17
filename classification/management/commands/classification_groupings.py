from django.core.management import BaseCommand

from classification.models import Classification, ClassificationModification, ClassificationSummaryCalculator
from classification.models.classification_grouping import ClassificationGrouping, AlleleOriginGrouping


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--summary', required=False, action="store_true")
        parser.add_argument('--all', required=False, action="store_true")
        parser.add_argument('--dirty', required=False, action="store_true")

    def handle(self, *args, **options):
        if not any(options.get(x) for x in ["summary", "all", "dirty"]):
            raise ValueError("Must provide one or more of summary, all, dirty")

        if options.get("summary"):
            for cm in ClassificationModification.objects.filter(is_last_published=True).select_related("classification"):
                classification = cm.classification
                classification.summary = ClassificationSummaryCalculator(cm).cache_dict()
                classification.save(update_fields=["summary"])

        if options.get("all"):
            for classification in Classification.objects.iterator():
                ClassificationGrouping.assign_grouping_for_classification(classification)
            ClassificationGrouping.objects.all().update(dirty=True)

        if options.get("all") or options.get("dirty"):
            for index, dirty in enumerate(ClassificationGrouping.objects.filter(dirty=True).iterator()):
                dirty.update()
                if index % 1000 == 0 and index:
                    print(f"Updating {index} classification groupings")

            for index, dirty in enumerate(AlleleOriginGrouping.objects.filter(dirty=True).iterator()):
                dirty.update()
                if index % 1000 == 0 and index:
                    print(f"Updating {index} allele groupings")
