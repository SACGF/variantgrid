"""
Deletes the sample/cohort annotated CSV exports written as .csv.zip, so the VCF and sample pages offer to
generate the plain .csv the export now writes (#2368). The pre_delete signal unlinks each file.
"""
import logging

from django.core.management.base import BaseCommand

from snpdb.models import CachedGeneratedFile

ANNOTATED_EXPORT_GENERATORS = ["export_cohort_to_downloadable_file", "export_sample_to_downloadable_file"]


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        zipped_qs = CachedGeneratedFile.objects.filter(generator__in=ANNOTATED_EXPORT_GENERATORS,
                                                       filename__endswith=".csv.zip")
        num_deleted, _ = zipped_qs.delete()
        logging.info("Deleted %d zipped annotated CSV exports", num_deleted)
