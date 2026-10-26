import logging
import os

from django.conf import settings

from annotation.annotation_versions import vav_diff_vs_kwargs
from annotation.celery_utils import annotation_worker_slots
from annotation.models import VariantAnnotationVersion
from annotation.vep_annotation import get_vep_variant_annotation_version_kwargs
from annotation.vep_plugin_backfill import vep_mismatch_fix
from snpdb.models import GenomeBuild

INSTALL_VEP = "https://github.com/SACGF/variantgrid/wiki/Install-VEP"


def vep_pool_oversubscription(worker_slots: int, fork: int, cpu_count: int) -> str | None:
    """ #2477: annotation_workers concurrency x ANNOTATION_VEP_FORK is how many CPU-bound VEPs the dispatcher runs
        at once when there is a backlog; more than the box has cores slows every run (and the web server) rather
        than finishing any sooner. None when the pool fits. """
    vep_processes = worker_slots * max(fork or 1, 1)
    if vep_processes <= cpu_count:
        return None
    return (f"annotation_workers concurrency ({worker_slots}) x ANNOTATION_VEP_FORK ({fork}) = {vep_processes} "
            f"VEP processes for {cpu_count} cores - lower --concurrency in config/celery/celeryd_annotation_workers.env "
            f"or ANNOTATION_VEP_FORK")


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
        vav = VariantAnnotationVersion.latest(genome_build)
        create_vav_cmd = f"python3 manage.py create_new_variant_annotation_version --genome-build={genome_build}"
        fix_vav = create_vav_cmd
        if vav is None:
            if vav_new := VariantAnnotationVersion.latest(genome_build,
                                                         status=VariantAnnotationVersion.Status.NEW):
                fix_vav = f"Your systems admin needs to promote VAV.pk={vav_new.pk} to ACTIVE."

        vep_data[f"have_variant_annotation_version_{genome_build}"] = {
            "valid": vav is not None,
            "fix": fix_vav,
        }

        run_vep_ok = False
        try:
            if vav.status == VariantAnnotationVersion.Status.ACTIVE:
                vep_vav_kwargs = get_vep_variant_annotation_version_kwargs(genome_build)
                run_vep_ok = True
                if vav:
                    # we only bother reporting diff vs latest if we can even run VEP (as they have bigger problems then)
                    diff = vav_diff_vs_kwargs(vav, vep_vav_kwargs)
                    warning = None
                    if diff:
                        warning = "Latest VariantAnnotationVersion doesn't match current VEP settings. " + \
                                  f"{diff=} {vep_mismatch_fix(vav, vep_vav_kwargs)}"
                    vep_data[f"latest_variant_annotation_{genome_build}_matches_current_vep"] = {
                        "valid": True,
                        "warning": warning,
                    }
        except Exception:
            pass

        logging.disable(previous_level)

        vep_data[f"run_vep_{genome_build}"] = {
            "valid": run_vep_ok,
            "fix": f"Try 'vep_version --genome-build={genome_build}' management command {INSTALL_VEP}",
        }

    # A warning, not invalid: an oversubscribed pool still finishes, slowly
    pool_check = {"valid": True}
    if warning := vep_pool_oversubscription(annotation_worker_slots(), settings.ANNOTATION_VEP_FORK,
                                            os.cpu_count() or 1):
        pool_check["warning"] = warning
    vep_data["vep_pool_fits_cores"] = pool_check

    return vep_data
