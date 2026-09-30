"""
Clears failed CachedGeneratedFile rows so the next request regenerates them. A failed generation stays
failed until the cause is fixed (@see snpdb.models.CachedGeneratedFile.needs_regenerating); the fix ships
with a migration that runs this command as a ManualOperation, filtered to the rows that failure produced.
"""
import logging

from django.core.management.base import BaseCommand, CommandError

from snpdb.models import CachedGeneratedFile


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument("--generator", help="eg export_cohort_to_downloadable_file")
        parser.add_argument("--exception-contains", help="Only rows whose stored exception contains this text")
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        generator = options["generator"]
        exception_contains = options["exception_contains"]
        if not (generator or exception_contains):
            raise CommandError("Provide --generator and/or --exception-contains - a fix clears the failures it "
                               "explains, not every failure")

        failed_qs = CachedGeneratedFile.objects.filter(exception__isnull=False)
        if generator:
            failed_qs = failed_qs.filter(generator=generator)
        if exception_contains:
            failed_qs = failed_qs.filter(exception__contains=exception_contains)

        for cgf in failed_qs:
            logging.info("%s %s: %s", "Would clear" if options["dry_run"] else "Clearing", cgf, cgf.exception)
        if not options["dry_run"]:
            num_deleted, _ = failed_qs.delete()
            logging.info("Cleared %d failed CachedGeneratedFile rows", num_deleted)
