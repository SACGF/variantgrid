# Can be removed once all systems migrated to new code written 11 Feb 2021
from django.core.management.base import BaseCommand

from library.guardian_utils import admin_bot
from library.utils import iter_fixed_chunks
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Variant, Allele


class Command(BaseCommand):
    def handle(self, *args, **options):
        for genome_build in GenomeBuild.builds_with_annotation():
            print(f"Handling {genome_build}")
            variant_qs = Variant.objects.filter(Variant.get_contigs_q(genome_build), varianttag__isnull=False)
            # Do in small chunks so we can save as we go - this already uses smaller batches internally
            chunk_size_kb = 2
            for variant_chunk in iter_fixed_chunks(variant_qs, chunk_size_kb * 1000):
                print(f"Handling {chunk_size_kb}k chunk")
                populate_clingen_alleles_for_variants(genome_build, variant_chunk)  # Will add VariantAlleles

            alleles = Allele.objects.filter(variantallele__variant__in=variant_qs)
            create_liftover_pipelines(admin_bot(), alleles, ImportSource.COMMAND_LINE, genome_build)
