"""
`vg inspect annotationrun <pk>`: one AnnotationRun - its VariantAnnotationVersion, pipeline type and
version, status, the range lock (variant pk range), counts at each stage, timings, the lease and task,
the error, and whether `vg status` counts it as in flight or as abandoned (nothing will dispatch it).
"""
from typing import Any

from django.utils import timezone

from annotation.abandoned_annotation_runs import abandoned_annotation_runs_by_reason
from annotation.enums import AnnotationStatus
from annotation.models.models_annotation_run import AnnotationRun
from library.vg.inspect import ref

TEXT_CHARS = 400
STDERR_LINES = 5
_TIMINGS = ("created", "dump_start", "dump_end", "annotation_start", "annotation_end", "upload_start", "upload_end")


def load(key: str) -> AnnotationRun:
    try:
        return AnnotationRun.objects.select_related("annotation_range_lock__version", "pipeline_version").get(pk=int(key))
    except (ValueError, AnnotationRun.DoesNotExist) as e:
        raise LookupError(f"No AnnotationRun with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    run = load(key)
    lock = run.annotation_range_lock
    incomplete = run.status not in AnnotationStatus.get_completed_states()
    abandoned = _abandoned_reason(run) if incomplete else None
    data: dict[str, Any] = {
        "id": run.pk,
        "status": run.get_status_display(),
        "in_flight": incomplete and abandoned is None,
        "abandoned": abandoned,
        "pipeline_type": run.get_pipeline_type_display(),
        "pipeline_version": str(run.pipeline_version) if run.pipeline_version_id else None,
        "external": run.external,
        "version": _version(lock.version) if lock else None,
        "range_lock": {"id": lock.pk, "min_variant": {"kind": "variant", "id": lock.min_variant_id},
                       "max_variant": {"kind": "variant", "id": lock.max_variant_id}, "count": lock.count} if lock else None,
        "counts": {"to_process": run.count, "dumped": run.dump_count, "annotated": run.annotated_count,
                   "vep_skipped": run.vep_skipped_count},
        "timings": {name: _when(value) for name in _TIMINGS if (value := getattr(run, name))},
        "task": {"task_id": run.task_id, "leased_by": run.leased_by, "lease_expires": _when(run.lease_expires),
                 "lease_live": bool(run.lease_expires and run.lease_expires >= timezone.now()),
                 "attempts": run.attempt_count, "dispatches": run.dispatch_count, "upload_attempts": run.upload_attempts},
        "error": _truncate(run.error_exception),
        "vep_warnings": _truncate(run.vep_warnings),
        "files": {"dump": run.vcf_dump_filename, "annotated": run.vcf_annotated_filename,
                  "vep_skipped": run.vep_skipped_variants_filename},
        "url": run.get_absolute_url(),
    }
    if depth >= 2:
        data["pipeline_command"] = _truncate(run.pipeline_command)
        data["stderr_tail"] = run.pipeline_stderr.strip().splitlines()[-STDERR_LINES:] if run.pipeline_stderr else []
    return data


def _version(vav) -> dict[str, Any]:
    return {**ref("annotationversion", vav), "build": vav.genome_build_id, "vep": vav.vep, "status": vav.get_status_display()}


def _abandoned_reason(run: AnnotationRun) -> str | None:
    """ The `vg status` abandoned bucket this run falls in, if any """
    runs = AnnotationRun.objects.filter(pk=run.pk)
    return next((reason for reason, qs in abandoned_annotation_runs_by_reason(runs).items() if qs.exists()), None)


def _when(value) -> str | None:
    return value.isoformat(timespec="seconds") if value else None


def _truncate(text: str | None) -> str | None:
    if not text:
        return None
    return text if len(text) <= TEXT_CHARS else text[:TEXT_CHARS] + f"… ({len(text):,} chars)"
