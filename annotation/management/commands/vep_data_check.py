import os
from itertools import pairwise

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.enums import VariantAnnotationPipelineType
from annotation.vep_annotation import VEPConfig, get_vep_command
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    category = "ops"

    def handle(self, *args, **options):

        for genome_build in GenomeBuild.builds_with_annotation():
            check_files = {}
            for key in ["reference_fasta", "cytoband"]:
                check_files[key] = genome_build.settings[key]

            vep_config = VEPConfig(genome_build)
            for key, rel_path in vep_config.vep_data.items():
                if rel_path is not None:
                    check_files[key] = vep_config[key]

            for key, filename in check_files.items():
                if not os.path.exists(filename):
                    print(f"{key}: {filename} MISSING")

            # A plugin can be newer than the VEP_plugins checkout (AVI.pm, #2426), which VEP only reports
            # once a run fails
            vep_command = get_vep_command("in.vcf", "out.vcf", genome_build, genome_build.annotation_consortium,
                                          VariantAnnotationPipelineType.STANDARD)
            for flag, plugin_arg in pairwise(vep_command):
                if flag == "--plugin":
                    plugin_name = plugin_arg.split(",")[0]
                    plugin_filename = os.path.join(settings.ANNOTATION_VEP_PLUGINS_DIR, f"{plugin_name}.pm")
                    if not os.path.exists(plugin_filename):
                        print(f"{genome_build} plugin {plugin_name}: {plugin_filename} MISSING")

        for data in settings.VCF_IMPORT_COMMON_FILTERS.values():
            key = "gnomad_af_filename"
            if filename := data.get(key):
                if not filename.startswith("/"):
                    filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, filename)
                    if not os.path.exists(filename):
                        print(f"{key}: {filename} MISSING")
