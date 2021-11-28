import logging

from django.core.management.base import BaseCommand

from snpdb.models import VCF, VariantZygosityCountForVCF, VariantZygosityCountForSample, VariantZygosityCountCollection
from snpdb.models.models_enums import ImportStatus
from snpdb.variant_zygosity_count import create_variant_zygosity_counts, update_variant_zygosity_count_for_vcf


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--name', required=True, help="Variant Zygosity Count Collection name")
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        name = options.get("name")
        collection = VariantZygosityCountCollection.objects.get(name=name)

        if options.get("clear"):
            existing_vcf_counts = VariantZygosityCountForVCF.objects.filter(collection=collection)
            existing_sample_counts = VariantZygosityCountForSample.objects.filter(collection=collection)

            if existing_vcf_counts.exists() or existing_sample_counts.exists():
                logging.info("Clearing existing global counts")
                VariantZygosityCountForVCF.objects.filter(collection=collection).delete()
                VariantZygosityCountForSample.objects.filter(collection=collection).delete()
                logging.info("Resetting to 0")
                collection.variantzygositycount_set.all().update(ref_count=0, het_count=0, hom_count=0)
                logging.info("Done")
            else:
                logging.info("No existing VCF/Sample counts exist - skipping clear")

        create_variant_zygosity_counts()
        for vcf in VCF.objects.filter(import_status=ImportStatus.SUCCESS):
            update_variant_zygosity_count_for_vcf(collection, vcf, "+")
