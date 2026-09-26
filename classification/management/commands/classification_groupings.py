from contextlib import nullcontext

from django.core.management import BaseCommand
from django.db import transaction

from classification.models import (
    Classification,
    ClassificationModification,
    ClassificationSummaryCalculator,
)
from classification.models.classification_grouping import (
    AlleleOriginGrouping,
    ClassificationGrouping,
)
from classification.services.overlaps_services import OverlapServices
from classification.signals.classification_hooks_grouping_search_terms import latest_annotation_versions_cached
from library.utils.collection_utils import batch_iterator


class Command(BaseCommand):
    category = "maintenance"

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

        # a full rebuild would otherwise notify labs of every existing discordance
        with OverlapServices.discordance_notifications_suppressed() if (all or refresh) else nullcontext():
            self._rebuild(summary=summary, all=all, dirty=dirty, refresh=refresh)

    @staticmethod
    def _rebuild(summary: bool, all: bool, dirty: bool, refresh: bool):
        if all or summary:
            modification_qs = ClassificationModification.objects.filter(is_last_published=True).select_related("classification")
            for batch in batch_iterator(modification_qs.iterator(), 1000):
                classifications = []
                for cm in batch:
                    classification = cm.classification
                    classification.summary = ClassificationSummaryCalculator(cm).cache_dict()
                    classifications.append(classification)
                Classification.objects.bulk_update(classifications, fields=["summary"])
                print(f"Updated {len(classifications)} classification summaries")

        if all:
            classification_qs = Classification.objects.select_related("lab", "allele_info__allele")
            for index, classification in enumerate(classification_qs.iterator()):
                # every grouping is marked dirty and updated below
                ClassificationGrouping.assign_grouping_for_classification(classification, force_dirty_up=False, update_new_grouping=False)
                if index % 1000 == 0 and index:
                    print(f"Updating {index} classification assigned to classification groups")
            print("About to update all classification groups")
            ClassificationGrouping.objects.all().update(dirty=True)

        if all or dirty or refresh:
            qs = ClassificationGrouping.objects.all()
            if not refresh:
                qs = qs.filter(dirty=True)

            grouping_ids = list(qs.order_by("pk").values_list("pk", flat=True))
            with OverlapServices.overlap_recalcs_deferred(), latest_annotation_versions_cached():
                for index, batch_ids in enumerate(batch_iterator(grouping_ids, 500)):
                    with transaction.atomic():
                        for grouping in ClassificationGrouping.objects.filter(pk__in=batch_ids):
                            grouping.update()
                    print(f"Updated {index * 500 + len(batch_ids)} of {len(grouping_ids)} classification groupings")
                print("Recalculating overlaps")
