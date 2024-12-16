import os

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from genes.models import TranscriptVersion
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
    active_builds = list(GenomeBuild.builds_with_annotation())
    active_build_names = set(str(gb) for gb in active_builds)

    for genome_build in active_builds:
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

        if settings.LIFTOVER_BCFTOOLS_ENABLED:
            annotation_build_config = settings.ANNOTATION[genome_build.name]
            for dest_genome_build, chain_filename in annotation_build_config["liftover"].items():
                if dest_genome_build in active_build_names:
                    key = f"bcftools_chain_{genome_build.name}_to_{dest_genome_build}"
                    check_files[key] = chain_filename

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

    return annotation_data


def check_cdot_data() -> dict:
    cdot_checks = {}

    for genome_build in GenomeBuild.builds_with_annotation():
        cdot_checks[f"cdot_{genome_build}"] = {
            "valid": TranscriptVersion.objects.filter(genome_build=genome_build).exists(),
            "fix": f"python manage import_gene_annotation --genome-build={genome_build.name}",
        }

    try:
        # Check that latest exists
        from cdot.data_release import get_latest_data_release_tag_name, _get_version_from_tag_name

        tag_name = get_latest_data_release_tag_name()
        cdot_data_version = _get_version_from_tag_name(tag_name, data_version=True)
        valid = False
        if last_tv := TranscriptVersion.objects.all().order_by("pk").last():
            our_latest_cdot = last_tv.data.get("cdot")
            valid = cdot_data_version == our_latest_cdot
        cdot_data = {
            "valid": valid,
            "notes": f"data version = latest ({cdot_data_version})",
            "fix": "python3 manage.py import_cdot_latest"
        }
        cdot_checks["latest_cdot_data"] = cdot_data
    except ImportError:
        # Will already be covered in library version > 0.2.26
        pass

    return cdot_checks