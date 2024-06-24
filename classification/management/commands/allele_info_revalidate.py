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

        genome_build_matchers: dict[GenomeBuild, HGVSMatcher] = {}
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            genome_build_matchers[genome_build] = HGVSMatcher(genome_build)

        for row_index, allele_info in enumerate(ImportedAlleleInfo.objects.iterator()):
            if row_index % 1000 == 0:
                print(f"Processed {row_index} rows")
            try:
                use_hgvs = allele_info.imported_c_hgvs or allele_info.imported_g_hgvs
                hgvs_matcher = genome_build_matchers[allele_info.imported_genome_build_patch_version.genome_build]
                vc_extra = hgvs_matcher.get_variant_coordinate_used_transcript_kind_method_and_matches_reference(use_hgvs)
                variant_coordinate = vc_extra.variant_coordinate

                if str(variant_coordinate) != allele_info.variant_coordinate:
                    allele_info.dirty_message = f"New variant coordinate: {str(variant_coordinate)}"
                    allele_info.save()
                else:
                    pass
            except Exception as ex:
                if allele_info.variant_coordinate:
                    allele_info.dirty_message = f"Could not re-resolve to variant coordinate despite currently being resolved: {ex}"
                    allele_info.save()
        print(f"Processed {row_index} rows")
