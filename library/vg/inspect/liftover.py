"""
`vg inspect liftover <pk>`: a LiftoverRun - source and destination build, conversion tool, who and when,
its allele source and the upload pipeline that ran it, AlleleLiftover counts by status, and the first few
failed and other allele liftovers (status, error) as allele refs.
"""
from typing import Any

from django.db.models import Count

from library.vg.inspect import capped, ref
from library.vg.inspect.common import allele_summary
from snpdb.enums import ProcessingStatus
from snpdb.models.models_variant import AlleleLiftover, LiftoverRun
from upload.models.models import UploadPipeline
from upload.models.models_uploaded_files import UploadedLiftover


def load(key: str) -> LiftoverRun:
    try:
        return LiftoverRun.objects.select_related("user").get(pk=int(key))
    except (ValueError, LiftoverRun.DoesNotExist) as e:
        raise LookupError(f"No LiftoverRun with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    run = load(key)
    allele_liftovers = AlleleLiftover.objects.filter(liftover=run)
    by_status = allele_liftovers.values_list("status").annotate(n=Count("pk")).order_by("status")
    data: dict[str, Any] = {
        "id": run.pk,
        "name": str(run),
        "source_build": run.source_genome_build_id,
        "build": run.genome_build_id,
        "tool": run.get_conversion_tool_display(),
        "user": run.user.username,
        "created": run.created.isoformat(timespec="seconds"),
        "allele_source": _allele_source(run),
        "source_vcf": run.source_vcf,
        "upload": _upload(run),
        "allele_liftovers": {ProcessingStatus(status).label: n for status, n in by_status},
        "url": run.get_absolute_url(),
    }
    if depth >= 2:
        rows = allele_liftovers.select_related("allele__clingen_allele")
        data["failed"] = capped(rows.filter(status=ProcessingStatus.ERROR).order_by("pk"), _allele_liftover)
        data["other_alleles"] = capped(rows.exclude(status=ProcessingStatus.ERROR).order_by("pk"), _allele_liftover)
    return data


def _allele_source(run: LiftoverRun) -> str | None:
    if run.allele_source_id is None:
        return None
    source = run.get_allele_source()
    return f"{type(source).__name__} {source.pk}"


def _upload(run: LiftoverRun) -> dict[str, Any] | None:
    uploaded = UploadedLiftover.objects.filter(liftover=run).first()
    pipeline = UploadPipeline.objects.filter(file_upload_id=uploaded.file_upload_id).first() if uploaded else None
    return ref("upload", pipeline, pipeline.get_status_display()) if pipeline else None


def _allele_liftover(allele_liftover: AlleleLiftover) -> dict[str, Any]:
    info = {**allele_summary(allele_liftover.allele), "status": allele_liftover.get_status_display()}
    if allele_liftover.error:
        error = allele_liftover.error.get(AlleleLiftover.ERROR_JSON_MESSAGE_KEY, allele_liftover.error)
        info["error"] = str(error)[:160]
    return info
