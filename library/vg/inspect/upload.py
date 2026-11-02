"""
`vg inspect upload <pk>`: an UploadPipeline - the file that went in (name, path, user, type), the pipeline's
status, progress and wall time, its top-level steps in run order (status, items, seconds, error) and what
it produced (the UploadData row and the object it points at: a VCF, a liftover run, a gene list...).
"""
from typing import Any

from library.vg.inspect import capped, ref
from snpdb.models.models_variant import LiftoverRun
from snpdb.models.models_vcf import VCF
from upload.models.models import UploadPipeline, UploadStep
from upload.uploaded_file_type import get_upload_data_for_uploaded_file

STEP_CAP = 30
ERROR_CHARS = 160
# Produced objects with an inspector of their own, keyed by model
_PRODUCED_KINDS = {VCF: "vcf", LiftoverRun: "liftover"}


def load(key: str) -> UploadPipeline:
    try:
        return UploadPipeline.objects.select_related("file_upload__user").get(pk=int(key))
    except (ValueError, UploadPipeline.DoesNotExist) as e:
        raise LookupError(f"No UploadPipeline with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    pipeline = load(key)
    file_upload = pipeline.file_upload
    data: dict[str, Any] = {
        "id": pipeline.pk,
        "file": {"id": file_upload.pk, "name": file_upload.name, "type": file_upload.get_file_type_display(),
                 "source": file_upload.get_import_source_display(), "path": file_upload.get_filename(),
                 "exists": file_upload.exists},
        "user": file_upload.user.username,
        "created": file_upload.created.isoformat(timespec="seconds"),
        "status": pipeline.get_status_display(),
        "progress": f"{_progress_status(pipeline.progress_status)} ({pipeline.progress_percent:.0f}%)",
        "items_processed": pipeline.items_processed,
        "wall_seconds": pipeline.processing_seconds_wall_time,
        "cpu_seconds": pipeline.processing_seconds_cpu_time,
        "celery_task": pipeline.celery_task,
        "error_summary": pipeline.error_summary,
        "url": pipeline.get_absolute_url(),
        "produced": _produced(file_upload),
    }
    steps = UploadStep.objects.filter(upload_pipeline=pipeline)
    data["errors"] = [f"{name}: {error[:ERROR_CHARS]}"
                      for name, error in steps.exclude(error_message="").order_by("-pk").values_list("name", "error_message")[:3]]
    if depth >= 2:
        top_level = steps.filter(parent_upload_step__isnull=True).order_by("sort_order", "pk")
        data["steps"] = capped(top_level, _step, cap=STEP_CAP)
    else:
        data["steps"] = steps.count()
    return data


def _produced(file_upload) -> dict[str, Any] | None:
    upload_data = get_upload_data_for_uploaded_file(file_upload)
    if upload_data is None:
        return None
    produced: dict[str, Any] = {"upload_data": f"{type(upload_data).__name__} {upload_data.pk}",
                                "requires_user_input": upload_data.requires_user_input}
    if (obj := upload_data.get_data()) is not None:
        produced["data"] = ref(_PRODUCED_KINDS.get(type(obj), obj._meta.model_name), obj)
    return produced


def _progress_status(progress_status: str | None) -> str:
    """ A failed pipeline stores its traceback here: the last line is the exception """
    if not progress_status:
        return "-"
    lines = progress_status.strip().splitlines()
    return lines[0] if len(lines) == 1 else f"{lines[0]} ... {lines[-1][:ERROR_CHARS]}"


def _step(step: UploadStep) -> dict[str, Any]:
    seconds = (step.end_date - step.start_date).total_seconds() if step.start_date and step.end_date else None
    info: dict[str, Any] = {"name": step.name, "status": step.get_status_display(), "items": step.items_processed,
                            "seconds": round(seconds) if seconds is not None else None}
    if substeps := step.substep_set.count():
        info["substeps"] = substeps
    if step.error_message:
        info["error"] = step.error_message[:ERROR_CHARS]
    return info
