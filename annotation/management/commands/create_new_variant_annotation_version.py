import logging

from django.core.management.base import BaseCommand

from annotation.annotation_versions import get_or_create_variant_annotation_version_from_current_vep
from annotation.models import (
    AnnotationVersion,
    GeneAnnotationVersion,
    InvalidAnnotationVersionError,
    VariantAnnotationVersion,
)
from snpdb.models.models_genome import GenomeBuild

INSTALL_RELEASE_COMMAND = "python3 manage.py import_cdot_gene_annotation_release --genome-build="


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

            self._report_gene_annotation_release(vav)

    @staticmethod
    def _report_gene_annotation_release(vav: VariantAnnotationVersion):
        """ A build's gene set often doesn't change between VEP versions (GRCh37 in particular), so an
            existing release usually still matches - link that rather than making them install it again """
        if vav.gene_annotation_release is None:
            if not vav.link_gene_annotation_release():
                release_token = vav.cdot_gene_release_token or "could not be determined from this VEP"
                print(f"{vav.genome_build}: no GeneAnnotationRelease for "
                      f"{vav.get_annotation_consortium_display()} release '{release_token}'. "
                      f"To download and install it, run:")
                print(f"    {INSTALL_RELEASE_COMMAND}{vav.genome_build}")
                return
            print(f"{vav.genome_build}: linked existing GeneAnnotationRelease '{vav.gene_annotation_release}'")

        # The AnnotationVersion was built before we linked the release, so it can be pointing at a
        # GeneAnnotationVersion for a different release. Say so rather than leaving it to fail later on a
        # variant page.
        release = vav.gene_annotation_release
        if not GeneAnnotationVersion.objects.filter(gene_annotation_release=release).exists():
            print(f"{vav.genome_build}: GeneAnnotationRelease '{release}' has no GeneAnnotationVersion yet. "
                  f"To create it and bring the AnnotationVersion into line, run:")
            print(f"    {INSTALL_RELEASE_COMMAND}{vav.genome_build}")
            return

        av = AnnotationVersion.latest(vav.genome_build, validate=False, status=vav.status)
        if av is None:
            return
        try:
            av.validate()
        except InvalidAnnotationVersionError as e:
            print(f"{vav.genome_build}: AnnotationVersion {av.pk} is inconsistent - {e}")
            print(f"    {INSTALL_RELEASE_COMMAND}{vav.genome_build}")
