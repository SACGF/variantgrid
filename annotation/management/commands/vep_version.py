from django.core.management.base import BaseCommand

from annotation.vep_annotation import get_vep_version, VEPConfig, vep_dict_to_variant_annotation_version_kwargs
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', required=True)

    def handle(self, *args, **options):
        build_name = options["genome_build"]
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        vep_config = VEPConfig(genome_build)
        vep_version = get_vep_version(genome_build, vep_config.annotation_consortium)
        print("*" * 40)
        print("VEP kwargs:")
        print(vep_version)

        vav_kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_version)
        print("*" * 40)
        print("VariantAnnotationVersion kwargs:")
        print(vav_kwargs)
