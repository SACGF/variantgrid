"""
A TSO 500 CombinedVariantOutput used to import as a VCF of no records, whose one sample existed only
so the pair's patient chain had something to hang off. Nothing needs that Sample any more (#1904), so
the import writes no VCF - this soft-deletes the ones existing deployments already have, the way a
user deleting a VCF does, so they leave the data page and the run page.

The CombinedVariantOutput rows' rna_sample is SET_NULL, and patients.tasks.extraction_matching_tasks
refills it from the arm's real VCF.
"""
import logging

from django.core.management import BaseCommand

from snpdb.models import VCF
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs

SOURCE_PREFIX = "DRAGEN TSO500 CombinedVariantOutput"


def record_less_combined_variant_output_vcfs_qs():
    return VCF.objects.filter(source__startswith=SOURCE_PREFIX)


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        for vcf in record_less_combined_variant_output_vcfs_qs():
            logging.info("Soft deleting record-less CombinedVariantOutput VCF %s", vcf)
            soft_delete_vcfs(vcf.user, vcf.pk)
