from django.core.management.base import BaseCommand

from classification.models import ReclassificationEvent, ReclassificationEventBuilder


class Command(BaseCommand):
    """ Builds the clinical significance timelines that drive /classification/reclassification_analytics """

    def add_arguments(self, parser):
        parser.add_argument('--rebuild', action='store_true',
                            help="Recompute every classification's timeline, rather than only those without one")

    def handle(self, *args, **options):
        classification_qs = ReclassificationEventBuilder.tracked_classifications_qs()
        if not options["rebuild"]:
            classification_qs = classification_qs.exclude(
                pk__in=ReclassificationEvent.objects.values('classification_id'))

        classification_count = classification_qs.count()
        self.stdout.write(f"Building timelines for {classification_count} classifications")

        rows_written = ReclassificationEventBuilder.rebuild(classification_qs, progress=self._progress(classification_count))

        self.stdout.write(f"Wrote {rows_written} events, "
                          f"the timeline now holds {ReclassificationEvent.reclassifications_qs().count()} reclassifications")

    def _progress(self, classification_count: int):
        processed = 0
        while True:
            processed += 1
            if processed % 5000 == 0:
                self.stdout.write(f"... {processed} / {classification_count}")
            yield processed
