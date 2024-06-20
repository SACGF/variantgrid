from django.core.management import BaseCommand

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSMatcher
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        ###
        # Right now, just validated ImportedAlleleInfo variant-coordinate
        # In future will handle variant resolution too
        ###

        genome_build_matchers: dict[GenomeBuild, HGVSMatcher] = {}
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            genome_build_matchers[genome_build] = HGVSMatcher(genome_build)

        for allele_info in ImportedAlleleInfo.objects.all():
            try:
                use_hgvs = allele_info.imported_c_hgvs or allele_info.imported_g_hgvs
                hgvs_matcher = genome_build_matchers[allele_info.imported_genome_build_patch_version.genome_build]
                vc_extra = hgvs_matcher.get_variant_coordinate_used_transcript_kind_method_and_matches_reference(use_hgvs)
                variant_coordinate = vc_extra.variant_coordinate

                if str(variant_coordinate) != allele_info.variant_coordinate:
                    print(f"Allele Info {allele_info.pk} resolved to different variant coordinate. Old = {allele_info.variant_coordinate} New = {str(variant_coordinate)}")
            except:
                if allele_info.variant_coordinate:
                    print(f"Allele Info {allele_info.pk} could not resolve, though it currently has a cached variant coordinate")
