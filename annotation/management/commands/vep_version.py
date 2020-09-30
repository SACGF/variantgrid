from django.core.management.base import BaseCommand

from annotation.vep_annotation import get_vep_version, VEPConfig
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', required=True)

    def handle(self, *args, **options):
        build_name = options["genome_build"]
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        vep_config = VEPConfig(genome_build)
        vep_version = get_vep_version(genome_build, vep_config.annotation_consortium)
        print(vep_version)
