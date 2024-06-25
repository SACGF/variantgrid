import itertools

from django.core.management import BaseCommand

from classification.models import ImportedAlleleInfo
from classification.models.variant_resolver import VariantResolver
from genes.hgvs import HGVSMatcher
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        self.iterate_alleles()


    def iterate_alleles(self):
        ###
        # Right now, just validated ImportedAlleleInfo variant-coordinate
        # In future will handle variant resolution too
        ###
        for row_index, allele_info in enumerate(ImportedAlleleInfo.objects.filter(imported_c_hgvs__isnull=False).iterator()):
            if row_index % 1000 == 0:
                print(f"Processed {row_index} rows")

            allele_info.dirty_check()

        print(f"Processed {row_index} rows")
