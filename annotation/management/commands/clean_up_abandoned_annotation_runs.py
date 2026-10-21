"""
Delete or ERROR the AnnotationRuns nothing will ever dispatch - no range lock, or on a HISTORICAL
VariantAnnotationVersion (#1940). See annotation/abandoned_annotation_runs.py.
"""
from django.core.management.base import BaseCommand

from annotation.abandoned_annotation_runs import clean_up_abandoned_annotation_runs


class Command(BaseCommand):
    category = "ops"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would be deleted or marked error, changing nothing")

    def handle(self, *args, **options):
        outcomes = clean_up_abandoned_annotation_runs(dry_run=options["dry_run"])
        if not outcomes:
            self.stdout.write("No abandoned annotation runs")
        for outcome, count in sorted(outcomes.items()):
            self.stdout.write(f"{outcome}: {count}")
