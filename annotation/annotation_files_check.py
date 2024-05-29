import os

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from snpdb.models import GenomeBuild


def annotation_data_exists(flat=False):
    build_data = {}

    for genome_build in GenomeBuild.builds_with_annotation():
        check_files = {}
        file_exists = {}
        for key in ["reference_fasta", "cytoband"]:
            check_files[key] = genome_build.settings[key]

        if data := settings.VCF_IMPORT_COMMON_FILTERS.get(genome_build.name):
            key = "gnomad_af_filename"
            if filename := data.get(key):
                if not filename.startswith("/"):
                    filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, filename)
                check_files[key] = filename

        vep_config = VEPConfig(genome_build)
        for key, rel_path in vep_config.vep_data.items():
            if rel_path is not None:
                check_files[key] = vep_config[key]

        for key, filename in check_files.items():
            file_exists[key] = os.path.exists(filename)

        build_data[genome_build.name] = file_exists

    if flat:
        flat_data = {}
        for genome_build_name, file_exists in build_data.items():
            for key, exists in file_exists.items():
                flat_data[f"{genome_build_name}/{key}"] = exists
        build_data = flat_data

    return build_data
