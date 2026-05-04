from django.db.models import Count

from annotation.models import AnnotationRangeLock, AnnotationRun, VariantAnnotationVersion, AnnotationVersion
from annotation.models.models_enums import AnnotationStatus
from snpdb.models import GenomeBuild


def is_variant_annotation_version_populated(vav: VariantAnnotationVersion) -> bool:
    """ True when the VAV's annotation tables are populated:
        - At least one AnnotationRangeLock exists.
        - All AnnotationRuns for this VAV are FINISHED (none queued/running/errored). """
    arl_qs = AnnotationRangeLock.objects.filter(version=vav)
    if not arl_qs.exists():
        return False
    runs_qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
    if not runs_qs.exists():
        return False
    unfinished = runs_qs.exclude(status=AnnotationStatus.FINISHED).count()
    return unfinished == 0


def get_variant_annotation_progress(vav: VariantAnnotationVersion) -> dict:
    """ Counts of annotated vs total range-lock variants for the VAV. """
    runs_qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
    finished = runs_qs.filter(status=AnnotationStatus.FINISHED).count()
    total = runs_qs.count()
    return {"finished": finished, "total": total}


def check_annotation_versions() -> dict:
    annotation_versions_checks = {}
    for genome_build in GenomeBuild.builds_with_annotation():
        build_av = {}
        try:
            AnnotationVersion.latest(genome_build, validate=True)
            build_av["valid"] = True
        except Exception as e:
            build_av["valid"] = False
            build_av["fix"] = f"{e}: See 'Annotation' web page for details"
        annotation_versions_checks[f"Annotation Version for {genome_build=}"] = build_av
    return annotation_versions_checks


def check_variant_annotation_runs_status() -> dict:
    # I am going to make this a warning for a while - before making it an error
    ARL_DUPE_ERROR = False
    # see https://github.com/SACGF/variantgrid_shariant/issues/177
    arl_qs = AnnotationRangeLock.objects.values('version_id', 'min_variant_id', 'max_variant_id')
    arl_qs = arl_qs.annotate(id_count=Count('id')).filter(id_count__gt=1)
    arl_message = "Remove duplicate AnnotationRangeLock objects"
    has_dupes = arl_qs.exists()
    if ARL_DUPE_ERROR:
        arl_data = {
            "valid": has_dupes,
            "fix": arl_message,
        }
    else:
        arl_data = {
            "valid": True,  # Just a warning
        }
        if has_dupes:
            arl_data["warning"] = arl_message + " (this is a just a warning but will become an error in the future)"
    annotation_status = {
        "AnnotationRangeLock duplicates": arl_data,
    }

    for genome_build in GenomeBuild.builds_with_annotation():
        vav = VariantAnnotationVersion.latest(genome_build)
        num_not_success = AnnotationRun.count_not_successful_runs_for_version(vav)
        anno_data = {
            "valid": True,  # Just a warning
        }
        if num_not_success:
            anno_data["warning"] = f"There are {num_not_success} annotation runs for {vav} with status error/incomplete"
        annotation_status[f"variant_annotation_latest_{genome_build}"] = anno_data

    return annotation_status
