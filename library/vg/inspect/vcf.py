"""
`vg inspect vcf <pk>`: an imported VCF - who and when, build, import status and its upload pipeline
(steps, error), the samples, the automatic cohort and its genotype collection, filters and header
size, archive state.
"""
from typing import Any

from library.vg.inspect import capped, ref
from library.vg.inspect.common import sample_summary
from snpdb.models import VCF, Cohort
from upload.models.models import UploadedVCF, UploadStep


def load(key: str) -> VCF:
    try:
        return VCF.objects.select_related("genome_build", "user", "project").get(pk=int(key))
    except (ValueError, VCF.DoesNotExist) as e:
        raise LookupError(f"No VCF with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    vcf = load(key)
    data: dict[str, Any] = {
        "id": vcf.pk,
        "name": vcf.name,
        "date": vcf.date.date().isoformat(),
        "user": vcf.user.username,
        "build": vcf.genome_build_id,
        "import_status": vcf.get_import_status_display(),
        "source": vcf.source or None,
        "project": vcf.project.name if vcf.project_id else None,
        "genotype_samples": vcf.genotype_samples,
        "header_lines": len(vcf.header.splitlines()) if vcf.header else 0,
        "filters": list(vcf.vcffilter_set.order_by("filter_id").values_list("filter_id", flat=True)),
        "zygosity_counted": vcf.variant_zygosity_count,
        "archived": vcf.data_archived_date.isoformat() if getattr(vcf, "data_archived_date", None) else None,
        "url": vcf.get_absolute_url(),
    }
    data["upload"] = _upload(vcf)
    cohort = Cohort.objects.filter(vcf=vcf).first()
    data["cohort"] = _cohort(cohort) if cohort else None
    if depth >= 2:
        data["samples"] = capped(vcf.sample_set.order_by("pk"), sample_summary)
    return data


def _upload(vcf: VCF) -> dict[str, Any] | None:
    uploaded = UploadedVCF.objects.filter(vcf=vcf).select_related("upload_pipeline__file_upload", "vcf_importer").first()
    if uploaded is None:
        return None
    pipeline = uploaded.upload_pipeline
    info: dict[str, Any] = {"file": uploaded.file_upload.name if uploaded.file_upload_id else None,
                            "path": uploaded.file_upload.path if uploaded.file_upload_id else None,
                            "importer": f"{uploaded.vcf_importer.name} {uploaded.vcf_importer.version}" if uploaded.vcf_importer_id else None}
    if pipeline:
        steps = UploadStep.objects.filter(upload_pipeline=pipeline)
        errors = steps.exclude(error_message="").order_by("-pk").values_list("name", "error_message")[:3]
        info.update({"pipeline": pipeline.pk, "status": pipeline.get_status_display(), "steps": steps.count(),
                     "processed": pipeline.items_processed, "wall_seconds": pipeline.processing_seconds_wall_time,
                     "errors": [f"{name}: {error[:160]}" for name, error in errors]})
    return info


def _cohort(cohort: Cohort) -> dict[str, Any]:
    info = {**ref("cohort", cohort, cohort.name), "samples": cohort.sample_count, "version": cohort.version,
            "import_status": cohort.get_import_status_display()}
    try:
        collection = cohort.cohort_genotype_collection
        info["genotype_collection"] = {"id": collection.pk, "table": collection.get_partition_table(),
                                       "common_split": collection.common_collection_id is not None}
    except Exception as e:  # pylint: disable=broad-exception-caught
        info["genotype_collection"] = f"none ({type(e).__name__})"
    return info
