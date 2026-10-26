import os

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.vep_annotation import get_missing_vep_files
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    category = "ops"

    def handle(self, *args, **options):

        for genome_build in GenomeBuild.builds_with_annotation():
            for key, filename in get_missing_vep_files(genome_build).items():
                print(f"{genome_build} {key}: {filename} MISSING")

        for data in settings.VCF_IMPORT_COMMON_FILTERS.values():
            key = "gnomad_af_filename"
            if filename := data.get(key):
                if not filename.startswith("/"):
                    filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, filename)
                    if not os.path.exists(filename):
                        print(f"{key}: {filename} MISSING")
