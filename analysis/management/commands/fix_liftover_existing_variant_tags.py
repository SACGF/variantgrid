# Can be removed once all systems migrated to new code written 11 Feb 2021
from django.core.management.base import BaseCommand

from library.guardian_utils import admin_bot
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import VariantAlleleCollectionSource, GenomeBuild, ImportSource, VariantAllele, \
    VariantAlleleCollectionRecord, Variant


class Command(BaseCommand):
    def handle(self, *args, **options):
        for genome_build in GenomeBuild.builds_with_annotation():
            variant_qs = Variant.objects.filter(Variant.get_contigs_q(genome_build), varianttag__isnull=False)
            populate_clingen_alleles_for_variants(genome_build, variant_qs)  # Will add VariantAlleles

            va_collection = VariantAlleleCollectionSource.objects.create(genome_build=genome_build)
            records = []
            for va in VariantAllele.objects.filter(variant__in=variant_qs):  # VariantAlleles added above
                records.append(VariantAlleleCollectionRecord(collection=va_collection, variant_allele=va))

            if records:
                VariantAlleleCollectionRecord.objects.bulk_create(records, batch_size=2000)
            create_liftover_pipelines(admin_bot(), va_collection, ImportSource.COMMAND_LINE, genome_build)
