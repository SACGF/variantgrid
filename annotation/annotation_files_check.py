import os

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from snpdb.models import GenomeBuild

_VARIANTGRID_DOWNLOAD_BASE_DIR = "http://variantgrid.com/download/annotation/VEP"


def _get_fix_instructions(filename) -> str:
    dirname = os.path.dirname(filename)
    vg_path = filename.replace(settings.ANNOTATION_VEP_BASE_DIR, _VARIANTGRID_DOWNLOAD_BASE_DIR)
    fix_instructions = f"cd {dirname};wget {vg_path}"
    return fix_instructions


def annotation_data_exists(flat=False) -> dict:
    all_build_data = {}
    # We can sometimes use files twice, only report once
    unique_filenames = set()

    for genome_build in GenomeBuild.builds_with_annotation():
        check_files = {}
        file_data = {}
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

        for base_key, base_filename in check_files.items():
            if base_filename in unique_filenames:
                continue
            else:
                unique_filenames.add(base_filename)

            files = {base_key: base_filename}
            if ".vcf" in base_filename:
                files[f"{base_key}_tbi"] = base_filename + ".tbi"

            for key, filename in files.items():
                file_data[key] = {
                    "valid": os.path.exists(filename),
                    "fix": _get_fix_instructions(filename),
                }
        all_build_data[genome_build.name] = file_data

    if flat:
        flat_data = {}
        for genome_build_name, build_data in all_build_data.items():
            for key, data in build_data.items():
                flat_data[f"{genome_build_name}/{key}"] = data
        annotation_data = flat_data
    else:
        annotation_data = all_build_data


    # VEP dirs
    VEP_DIRS = {
        "vep_cache": settings.ANNOTATION_VEP_CACHE_DIR,
        "vep_plugins": settings.ANNOTATION_VEP_PLUGINS_DIR,
    }
    for k, dir_name in VEP_DIRS.items():
        valid = os.path.exists(dir_name)
        annotation_data[k] = {
            "valid": valid,
            "fix": "https://github.com/SACGF/variantgrid/wiki/Install-VEP",
        }

    return annotation_data
