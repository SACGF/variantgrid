import logging
import os

from django.conf import settings

from annotation.annotation_versions import vav_diff_vs_kwargs
from annotation.models import VariantAnnotationVersion
from annotation.vep_annotation import get_vep_variant_annotation_version_kwargs
from snpdb.models import GenomeBuild


INSTALL_VEP = "https://github.com/SACGF/variantgrid/wiki/Install-VEP"

def check_vep() -> dict:
    vep_data = {}

    VEP_DIRS = {
        "vep_code": settings.ANNOTATION_VEP_CODE_DIR,
        "vep_cache_base_dir": settings.ANNOTATION_VEP_CACHE_DIR,
        "vep_plugins": settings.ANNOTATION_VEP_PLUGINS_DIR,
    }
    for k, dir_name in VEP_DIRS.items():
        valid = os.path.exists(dir_name)
        vep_data[k] = {
            "valid": valid,
            "fix": INSTALL_VEP,
        }

    # Cache and plugins already checked in "annotation_files_check"
    # VEP code dir exists

    for genome_build in GenomeBuild.builds_with_annotation():
        # Disable all the annoying VEP logging
        previous_level = logging.root.manager.disable
        logging.disable(logging.CRITICAL)
        vav = None
        try:
            vav = VariantAnnotationVersion.latest(genome_build)
        except:
            pass

        create_vav_cmd = f"python3 manage.py create_new_variant_annotation_version --genome-build={genome_build}"
        vep_data[f"have_variant_annotation_version_{genome_build}"] = {
            "valid": vav is not None,
            "fix": create_vav_cmd,
        }

        try:
            vep_vav_kwargs = get_vep_variant_annotation_version_kwargs(genome_build)
            run_vep_ok = True
            if vav:
                # we only bother reporting diff vs latest if we can even run VEP (as they have bigger problems then)
                # This is only a warning as may not be an actual error - may want to just leave as old version
                diff = vav_diff_vs_kwargs(vav, vep_vav_kwargs)
                warning = None
                if diff:
                    warning = "Lastest VariantAnnotationVersion doesn't match current VEP settings " + \
                              f"(warning only as this may be what you want). {diff=} Upgrade via '{create_vav_cmd}'"
                vep_data[f"latest_variant_annotation_{genome_build}_matches_current_vep"] = {
                    "valid": True,
                    "warning": warning,
                }
        except:
            run_vep_ok = False

        logging.disable(previous_level)

        vep_data[f"run_vep_{genome_build}"] = {
            "valid": run_vep_ok,
            "fix": f"Try 'vep_version --genome-build={genome_build}' management command {INSTALL_VEP}",
        }

    return vep_data