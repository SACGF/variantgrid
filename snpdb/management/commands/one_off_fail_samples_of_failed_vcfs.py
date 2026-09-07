"""
Samples used to be released by the stats task independently of their VCF, so a VCF whose pipeline
then failed could be left ImportStatus.ERROR with samples still marked success. Sample status now
only moves with the VCF (@see snpdb/import_status.py) - this brings the existing rows into line.
"""
import logging

from django.core.management import BaseCommand

from snpdb.models import ImportStatus, Sample


def samples_of_failed_vcfs_qs():
    """ Samples still open for use on a VCF that failed - deletion states are left alone """
    return Sample.objects.filter(vcf__import_status=ImportStatus.ERROR) \
        .exclude(import_status__in=[ImportStatus.ERROR, *ImportStatus.DELETION_STATES])


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        updated = samples_of_failed_vcfs_qs().update(import_status=ImportStatus.ERROR)
        logging.info("Set %d samples of failed VCFs to %s", updated, ImportStatus.ERROR)
