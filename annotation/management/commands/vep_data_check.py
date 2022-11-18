import os

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.vep_annotation import VEPConfig
from snpdb.models import GenomeBuild


class Command(BaseCommand):
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

        for build_name, data in settings.VCF_IMPORT_COMMON_FILTERS.items():
            key = "gnomad_af_filename"
            if filename := data.get(key):
                if not filename.startswith("/"):
                    filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, filename)
                    if not os.path.exists(filename):
                        print(f"{key}: {filename} MISSING")
