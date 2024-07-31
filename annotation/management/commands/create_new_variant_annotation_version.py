import logging

from django.core.management.base import BaseCommand

from annotation.annotation_versions import get_or_create_variant_annotation_version_from_current_vep
from annotation.vep_annotation import get_vep_version, VEPConfig, vep_dict_to_variant_annotation_version_kwargs
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
