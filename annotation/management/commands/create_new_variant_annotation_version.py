import logging

from django.core.management.base import BaseCommand

from annotation.annotation_versions import get_or_create_variant_annotation_version_from_current_vep
from annotation.models import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--genome-build')

    def handle(self, *args, **options):
        if build_name := options.get("genome_build"):
            genome_build = GenomeBuild.get_name_or_alias(build_name)
            genome_builds = [genome_build]
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        for genome_build in genome_builds:
            vav, created = get_or_create_variant_annotation_version_from_current_vep(genome_build)
            if created:
                logging.info("Created: %s", vav)
            else:
                logging.info("Existing matches current VEP: %s", vav)

            if vav.gene_annotation_release is None:
                self._report_gene_annotation_release(vav)

    @staticmethod
    def _report_gene_annotation_release(vav: VariantAnnotationVersion):
        """ A build's gene set often doesn't change between VEP versions (GRCh37 in particular), so an
            existing release usually still matches - link it rather than making them install it again """
        if release := vav.link_gene_annotation_release():
            print(f"{vav.genome_build}: linked existing GeneAnnotationRelease '{release}'")
            return

        release_token = vav.cdot_gene_release_token or "could not be determined from this VEP"
        print(f"{vav.genome_build}: no GeneAnnotationRelease for {vav.get_annotation_consortium_display()} "
              f"release '{release_token}'. To download and install it, run:")
        print(f"    python3 manage.py import_cdot_gene_annotation_release --genome-build={vav.genome_build}")
