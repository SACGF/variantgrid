"""
External annotation runs (#1568): dump/import helpers shared by the annotation_external command.

The heavy VEP step can be run off-VM and the resulting annotated VCFs re-imported, and reused between a
database and its own clone for identical annotation runs. See claude/plans/1568_external_annotation_runs_plan.md.
"""
import json
import logging
import os

from django.conf import settings
from django.db.models import Q
from django.utils import timezone

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks.annotate_variants import dump_variants
from annotation.vep_annotation import get_vep_variant_annotation_version_kwargs
from library.utils.file_utils import mk_path_for_file
from snpdb.models import Variant

# Bump if the dump metadata layout (build_dump_metadata) changes incompatibly (#1568)
DUMP_METADATA_SCHEMA_VERSION = 1


class VariantIdAlignmentError(ValueError):
    """ Raised when an externally-annotated VCF's recorded range endpoints do not match the local Variant
        coordinates at the matched run's range-lock min/max ids - i.e. it was produced against a different
        or diverged database (variant ids past the clone split). See #1568 §6a. """


def variant_annotation_version_identity(variant_annotation_version: VariantAnnotationVersion) -> dict:
    """ Version-identity used to match an externally-annotated VCF back to a local VariantAnnotationVersion
        (#1568). The field set is the canonical VEP/annotation identity defined by
        get_vep_variant_annotation_version_kwargs() (the same fields vav_diff_vs_kwargs / the ##VEP= header
        check compare) - we read those fields off the stored VAV rather than maintaining a parallel list.
        genome_build / annotation_consortium live at the top level of the dump metadata. """
    vep_kwargs = get_vep_variant_annotation_version_kwargs(variant_annotation_version.genome_build)
    identity = {field: getattr(variant_annotation_version, field)
                for field in vep_kwargs if field not in ("genome_build", "annotation_consortium")}
    gene_annotation_release = variant_annotation_version.gene_annotation_release
    identity["gene_annotation_release"] = {
        "pk": variant_annotation_version.gene_annotation_release_id,
        "version": gene_annotation_release.version if gene_annotation_release else None,
    }
    return identity


def build_dump_metadata(annotation_run: AnnotationRun) -> dict:
    """ Self-describing sidecar metadata for an external annotation dump (#1568).
        Drives import-time matching (version identity + range coordinate strings) and the §6a variant-id
        alignment check. annotation_run_pk is recorded for human/debug only - it is NOT a cross-DB handle
        (range locks/versions can be created in a different order on a clone). """
    annotation_range_lock = annotation_run.annotation_range_lock
    variant_annotation_version = annotation_run.variant_annotation_version
    return {
        "schema": DUMP_METADATA_SCHEMA_VERSION,
        "site_name": settings.SITE_NAME,
        "annotation_run_pk": annotation_run.pk,
        "pipeline_type": annotation_run.pipeline_type,
        "genome_build": variant_annotation_version.genome_build.name,
        "annotation_consortium": variant_annotation_version.annotation_consortium,
        "variant_annotation_version": {
            "pk": variant_annotation_version.pk,
            **variant_annotation_version_identity(variant_annotation_version),
        },
        "range": {
            "min_variant": str(annotation_range_lock.min_variant.coordinate),
            "max_variant": str(annotation_range_lock.max_variant.coordinate),
            "count": annotation_range_lock.count,
        },
        "batch": {
            "annotation_vep_batch_min": settings.ANNOTATION_VEP_BATCH_MIN,
            "annotation_vep_batch_max": settings.ANNOTATION_VEP_BATCH_MAX,
        },
        "dump_count": annotation_run.dump_count,
    }


def write_dump_metadata(annotation_run: AnnotationRun, dump_dir=None) -> str:
    """ Write the sidecar metadata JSON next to the dump VCF, returning its path (#1568). """
    meta_filename = annotation_run.get_dump_metadata_filename(dump_dir=dump_dir)
    mk_path_for_file(meta_filename)
    with open(meta_filename, "w") as f:
        json.dump(build_dump_metadata(annotation_run), f, indent=2)
    return meta_filename


def parse_dump_metadata(path) -> dict:
    """ Read + validate a sidecar metadata JSON written by write_dump_metadata (#1568). """
    with open(path) as f:
        meta = json.load(f)
    schema = meta.get("schema")
    if schema != DUMP_METADATA_SCHEMA_VERSION:
        raise ValueError(f"{path}: dump metadata schema {schema!r} != expected {DUMP_METADATA_SCHEMA_VERSION}")
    return meta


def dump_external_annotation_runs(variant_annotation_version: VariantAnnotationVersion,
                                  output_dir: str,
                                  pipeline_type=VariantAnnotationPipelineType.STANDARD) -> list[AnnotationRun]:
    """ Create all range locks + external AnnotationRuns for the (NEW) version up front, dumping each VCF +
        sidecar metadata into output_dir and parking it in EXTERNAL_DUMP_COMPLETED awaiting off-VM annotation
        (#1568 §4.2).

        Range locks are sized by the target's ANNOTATION_VEP_BATCH_MIN/MAX so boundaries line up with what a
        clone would produce (§3) - external parallelism comes from forks/concurrency across runs, never larger
        locks. """
    os.makedirs(output_dir, exist_ok=True)
    annotation_runs = []
    while True:
        range_lock, _unannotated_count = get_annotation_range_lock_and_unannotated_count(
            variant_annotation_version, settings.ANNOTATION_VEP_BATCH_MIN, settings.ANNOTATION_VEP_BATCH_MAX)
        if range_lock is None:
            break
        range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=range_lock,
                                                      pipeline_type=pipeline_type, external=True)
        dump_count = dump_variants(annotation_run, dump_dir=output_dir)
        meta_filename = write_dump_metadata(annotation_run, dump_dir=output_dir)
        logging.info("Dumped external AnnotationRun %s: %d variants -> %s (meta %s)",
                     annotation_run.pk, dump_count, annotation_run.vcf_dump_filename, meta_filename)
        annotation_runs.append(annotation_run)

    logging.info("Dumped %d external annotation run(s) for %s into %s",
                 len(annotation_runs), variant_annotation_version, output_dir)
    return annotation_runs


def dump_existing_annotation_runs(variant_annotation_version: VariantAnnotationVersion,
                                  output_dir: str,
                                  pipeline_type=VariantAnnotationPipelineType.STANDARD,
                                  leave: int = 0) -> list[AnnotationRun]:
    """ Adopt AnnotationRuns that already exist for `variant_annotation_version` in CREATED state (created by
        the normal scheduler but not yet annotated), marking them external and dumping each VCF + sidecar
        metadata into output_dir (#1568). Unlike dump_external_annotation_runs (which creates every run up
        front for a NEW version), this leaves the existing range locks untouched and only claims runs the
        dispatcher has not already leased.

        `leave` keeps that many of the lowest-min-variant CREATED runs on the normal in-VM pipeline so
        annotation can be parallelised: the local VM keeps chewing through the low end (moving the unannotated
        watermark) while the dumped remainder is annotated off-VM. """
    if leave < 0:
        raise ValueError(f"leave must be >= 0, got {leave}")

    os.makedirs(output_dir, exist_ok=True)
    now = timezone.now()
    # Mirror the dispatcher's dispatchable filter (annotation_scheduler_task._dispatchable_runs_qs) and its
    # lowest-min-variant-first order, so we adopt exactly the runs it would otherwise launch.
    dispatchable = AnnotationRun.objects.filter(
        annotation_range_lock__version=variant_annotation_version,
        pipeline_type=pipeline_type,
        external=False,
        task_id__isnull=True,
        status=AnnotationStatus.CREATED,
    ).filter(Q(lease_expires__isnull=True) | Q(lease_expires__lt=now)) \
        .order_by("annotation_range_lock__min_variant_id")

    candidate_ids = list(dispatchable.values_list("pk", flat=True))
    kept, candidate_ids = candidate_ids[:leave], candidate_ids[leave:]
    logging.info("dump_existing: %d dispatchable run(s); leaving %d on the in-VM pipeline, dumping %d",
                 len(kept) + len(candidate_ids), len(kept), len(candidate_ids))

    annotation_runs = []
    for pk in candidate_ids:
        # Atomically claim as external only while still dispatchable (same filter as the dispatcher) so we
        # never adopt a run it just leased. If we lose the race (0 rows updated) skip it; if we win, the run
        # is external and annotate_variants no-ops even if the dispatcher launches it.
        claimed = AnnotationRun.objects.filter(
            pk=pk,
            external=False,
            task_id__isnull=True,
            status=AnnotationStatus.CREATED,
        ).filter(Q(lease_expires__isnull=True) | Q(lease_expires__lt=now)).update(external=True)
        if claimed != 1:
            logging.warning("Skipping AnnotationRun %s - no longer dispatchable (leased/launched by the "
                            "scheduler?)", pk)
            continue

        annotation_run = AnnotationRun.objects.get(pk=pk)
        dump_count = dump_variants(annotation_run, dump_dir=output_dir)
        meta_filename = write_dump_metadata(annotation_run, dump_dir=output_dir)
        logging.info("Dumped existing AnnotationRun %s: %d variants -> %s (meta %s)",
                     annotation_run.pk, dump_count, annotation_run.vcf_dump_filename, meta_filename)
        annotation_runs.append(annotation_run)

    logging.info("Dumped %d existing annotation run(s) for %s into %s",
                 len(annotation_runs), variant_annotation_version, output_dir)
    return annotation_runs


def verify_annotated_vcf_variant_ids(annotation_run: AnnotationRun, meta: dict):
    """ §6a variant-id alignment safety check. The variant_id baked into an external dump's INFO is the
        ORIGIN DB's Variant.id. On a pre-split clone those ids line up 1:1 with this DB, so the local Variant
        at the matched run's range-lock min/max id has the same coordinate the dump recorded. A post-split
        file (or one from the wrong DB) shifts at least one endpoint to a different/missing coordinate.

        Range locks are a contiguous block of Variant.pk and are deterministic in pk order + batch size (§3),
        so if both endpoints align the interior does too - checking endpoints is sufficient.

        Raises VariantIdAlignmentError with an actionable message on mismatch/missing variant; the import
        command catches this to fail only this run (not the whole import). """
    annotation_range_lock = annotation_run.annotation_range_lock
    meta_range = meta["range"]

    for endpoint, local_variant_id, expected_coordinate in [
        ("min", annotation_range_lock.min_variant_id, meta_range["min_variant"]),
        ("max", annotation_range_lock.max_variant_id, meta_range["max_variant"]),
    ]:
        local_variant = Variant.objects.filter(pk=local_variant_id).first()
        if local_variant is None:
            raise VariantIdAlignmentError(
                f"AnnotationRun {annotation_run.pk}: range {endpoint} local variant {local_variant_id} "
                f"does not exist but annotated file recorded {expected_coordinate!r} - produced against a "
                f"different/diverged database (id past the split?). Skipping this run; re-run it normally.")

        # Compute the coordinate string the same way the dump/meta wrote it (Variant.coordinate string repr)
        local_coordinate = str(local_variant.coordinate)
        if local_coordinate != expected_coordinate:
            raise VariantIdAlignmentError(
                f"AnnotationRun {annotation_run.pk}: range {endpoint} local variant {local_variant_id} is "
                f"{local_coordinate} but annotated file recorded {expected_coordinate} - produced against a "
                f"different/diverged database (id past the split?). Skipping this run; re-run it normally.")
