import os

import requests
from cdot.data_release import _get_version_from_tag_name, get_latest_data_release_tag_name
from django.conf import settings

from annotation.vep_annotation import (
    VEPConfig,
    get_vep_plugin_files,
    get_vep_plugin_url,
    is_variantgrid_vep_plugin,
)
from genes.models import CdotDataVersion, TranscriptVersion
from snpdb.models import GenomeBuild
from snpdb.models.models_somalier import SomalierConfig

_VARIANTGRID_DOWNLOAD_BASE_DIR = "https://variantgrid.com/download/annotation"


def _get_fix_instructions(filename) -> str:
    """'cd <dir>;wget <url>' for a file variantgrid.com hosts, anchored at ANNOTATION_BASE_DIR (not the narrower VEP
    subdir) so reference_fasta in fasta/ and the liftover chains map to a real URL. A file outside it ships with the
    repo (the cytobands under VARIANTGRID_REPO_REFERENCE_DIR), so the fix is a checkout, not a download."""
    if not filename.startswith(settings.ANNOTATION_BASE_DIR):
        return f"git checkout -- {os.path.relpath(filename, settings.BASE_DIR)}"
    dirname = os.path.dirname(filename)
    vg_path = filename.replace(settings.ANNOTATION_BASE_DIR, _VARIANTGRID_DOWNLOAD_BASE_DIR)
    return f"cd {dirname};wget {vg_path}"


def annotation_data_exists(flat=False, include_tbi_for_gz=False, include_csi=False) -> dict:
    """Per build (or flat as '<build>/<key>'), each annotation file as {"valid": exists, "fix": "cd ..;wget .."}.

    include_csi adds a '<key>_csi' entry, marked "optional": True, beside every .tbi: a fine-binned CSI index that
    htslib prefers over the .tbi and that makes VEP's custom/plugin lookups several times faster (#2499). It is
    optional because variantgrid.com only hosts one for the files it builds, and a deployment runs correctly on the
    .tbi alone; annotation/annotation_data/generate_annotation/build_tabix_csi.sh builds the rest locally."""
    all_build_data = {}
    # We can sometimes use files twice, only report once
    unique_filenames = set()
    active_builds = list(GenomeBuild.builds_with_annotation())
    active_build_names = {str(gb) for gb in active_builds}

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
            liftover_config = settings.ANNOTATION[genome_build.name]["liftover"]
            if liftover_fasta := liftover_config.get("fasta"):
                check_files["liftover_fasta"] = liftover_fasta
            for dest_genome_build, chain_filename in liftover_config["chain"].items():
                if dest_genome_build in active_build_names:
                    key = f"bcftools_chain_{genome_build.name}_to_{dest_genome_build}"
                    check_files[key] = chain_filename

        if settings.SOMALIER["enabled"]:
            # Not build specific - unique_filenames reports it under the first build only
            check_files["somalier"] = SomalierConfig().get_annotation("command")

        vep_config = VEPConfig(genome_build)
        for key, rel_path in vep_config.vep_data.items():
            if rel_path is not None:
                check_files[key] = vep_config[key]

        for base_key, base_filename in check_files.items():
            if base_filename in unique_filenames:
                continue

            unique_filenames.add(base_filename)
            files = {base_key: base_filename}
            vep_gz = (base_filename.endswith((".gz", ".bgz"))
                      and base_filename.startswith(settings.ANNOTATION_VEP_BASE_DIR))
            optional_keys = set()
            if ".vcf" in base_filename or (include_tbi_for_gz and vep_gz):
                files[f"{base_key}_tbi"] = base_filename + ".tbi"
            if include_csi and (".vcf" in base_filename or vep_gz):
                files[f"{base_key}_csi"] = base_filename + ".csi"
                optional_keys.add(f"{base_key}_csi")

            for key, filename in files.items():
                file_data[key] = {
                    "valid": os.path.exists(filename),
                    "fix": _get_fix_instructions(filename),
                }
                if key in optional_keys:
                    file_data[key]["optional"] = True

        # A plugin can be newer than the VEP_plugins checkout (AVI.pm, #2426), so comes from Ensembl rather than us.
        # Our own plugins are "refresh": the hosted copy is authoritative and replaces whatever is installed
        for plugin_name, filename in get_vep_plugin_files(genome_build).items():
            if filename in unique_filenames:
                continue
            unique_filenames.add(filename)
            file_data[f"plugin_{plugin_name}"] = {
                "valid": os.path.exists(filename),
                "fix": f"cd {os.path.dirname(filename)};wget {get_vep_plugin_url(plugin_name)}",
            }
            if is_variantgrid_vep_plugin(plugin_name):
                file_data[f"plugin_{plugin_name}"]["refresh"] = True
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
            "fix": f"python3 manage.py import_cdot_latest --genome-build={genome_build.name}",
        }

    # Only the build/consortium combinations a deployment has imported - one never imported is cdot_{build}'s job
    cdot_qs = CdotDataVersion.objects.filter(genome_build__in=GenomeBuild.builds_with_annotation())
    our_cdot_versions = set(cdot_qs.values_list("cdot_version", flat=True))
    try:
        tag_name = get_latest_data_release_tag_name()
    except requests.RequestException as e:
        # Unauthenticated GitHub API calls are limited to 60/hour per IP - not worth failing a deploy over
        cdot_checks["latest_cdot_data"] = {
            "valid": True,
            "warning": f"Could not check latest cdot data release (have {', '.join(sorted(our_cdot_versions))}): {e}",
        }
    else:
        cdot_data_version = _get_version_from_tag_name(tag_name, data_version=True)
        cdot_checks["latest_cdot_data"] = {
            "valid": our_cdot_versions == {cdot_data_version},
            "notes": f"data version ({', '.join(sorted(our_cdot_versions))}) = latest ({cdot_data_version})",
            "fix": "python3 manage.py import_cdot_latest"
        }

    return cdot_checks
