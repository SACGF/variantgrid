"""
External annotation runs (#1568): dump/import helpers shared by the annotation_external command.

The heavy VEP step can be run off-VM and the resulting annotated VCFs re-imported, and reused between a
database and its own clone for identical annotation runs. See claude/plans/1568_external_annotation_runs_plan.md.
"""
import glob
import json
import logging
import os
import shutil
import time
from typing import Optional

from django.conf import settings
from django.db.models import Q
from django.utils import timezone

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.models.models import AnnotationRangeLock, AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks.annotate_variants import dump_variants
from annotation.tasks.annotation_scheduler_task import dispatch_annotation_runs
from annotation.vep_annotation import get_vep_command, get_vep_variant_annotation_version_kwargs
from library.utils.file_utils import mk_path_for_file
from snpdb.models import GenomeBuild, Variant

# Bump if the dump metadata layout (build_dump_metadata) changes incompatibly (#1568)
DUMP_METADATA_SCHEMA_VERSION = 1

# Suffix the Snakemake bundle appends to each dump stem for the annotated VCF; import discovers annotated
# files alongside their .meta.json sidecar by this suffix (#1568).
ANNOTATED_VCF_SUFFIX = ".vep_annotated.vcf.gz"

# How often import_external_annotation_runs emits a running-progress line for large (e.g. 1k-run) imports.
IMPORT_PROGRESS_INTERVAL_SECONDS = 10


def _import_progress_line(processed: int, total: int, report: dict) -> str:
    tallies = " ".join(f"{key}={len(value)}" for key, value in report.items())
    return f"External annotation import: {processed}/{total} meta files ({tallies})"


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


def build_dump_metadata(annotation_run: AnnotationRun, identity: dict = None) -> dict:
    """ Self-describing sidecar metadata for an external annotation dump (#1568).
        Drives import-time matching (version identity + range coordinate strings) and the §6a variant-id
        alignment check. annotation_run_pk is recorded for human/debug only - it is NOT a cross-DB handle
        (range locks/versions can be created in a different order on a clone).

        `identity` is variant_annotation_version_identity() for this run's VAV; every run in a dump shares
        the same VAV so callers compute it once and pass it in (it runs VEP - see that function). """
    annotation_range_lock = annotation_run.annotation_range_lock
    variant_annotation_version = annotation_run.variant_annotation_version
    if identity is None:
        identity = variant_annotation_version_identity(variant_annotation_version)
    return {
        "schema": DUMP_METADATA_SCHEMA_VERSION,
        "site_name": settings.SITE_NAME,
        "annotation_run_pk": annotation_run.pk,
        "pipeline_type": annotation_run.pipeline_type,
        "genome_build": variant_annotation_version.genome_build.name,
        "annotation_consortium": variant_annotation_version.annotation_consortium,
        "variant_annotation_version": {
            "pk": variant_annotation_version.pk,
            **identity,
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


def write_dump_metadata(annotation_run: AnnotationRun, dump_dir=None, identity: dict = None) -> str:
    """ Write the sidecar metadata JSON next to the dump VCF, returning its path (#1568).
        Pass a precomputed `identity` (see build_dump_metadata) to avoid re-running VEP per run. """
    meta_filename = annotation_run.get_dump_metadata_filename(dump_dir=dump_dir)
    mk_path_for_file(meta_filename)
    with open(meta_filename, "w") as f:
        json.dump(build_dump_metadata(annotation_run, identity=identity), f, indent=2)
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
    # Every run in this dump shares the same VAV, so compute the version identity (which runs VEP) once,
    # lazily on the first run so a no-op dump never invokes VEP.
    identity = None
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
        if identity is None:
            identity = variant_annotation_version_identity(variant_annotation_version)
        meta_filename = write_dump_metadata(annotation_run, dump_dir=output_dir, identity=identity)
        logging.info("Dumped external AnnotationRun %s: %d variants -> %s (meta %s)",
                     annotation_run.pk, dump_count, annotation_run.vcf_dump_filename, meta_filename)
        annotation_runs.append(annotation_run)

    write_snakemake_bundle(output_dir, variant_annotation_version, pipeline_type=pipeline_type)
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

    # Every run in this dump shares the same VAV, so compute the version identity (which runs VEP) once,
    # lazily on the first claimed run so a no-op dump never invokes VEP.
    identity = None
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
        if identity is None:
            identity = variant_annotation_version_identity(variant_annotation_version)
        meta_filename = write_dump_metadata(annotation_run, dump_dir=output_dir, identity=identity)
        logging.info("Dumped existing AnnotationRun %s: %d variants -> %s (meta %s)",
                     annotation_run.pk, dump_count, annotation_run.vcf_dump_filename, meta_filename)
        annotation_runs.append(annotation_run)

    write_snakemake_bundle(output_dir, variant_annotation_version, pipeline_type=pipeline_type)
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


# --------------------------------------------------------------------------------------------------------
# Snakemake bundle generation (#1568 §4.3): emitted into the dump dir so the operator copies the whole
# directory to a compute box, edits config.yaml, and runs `snakemake` to produce annotated VCFs. Reusing
# get_vep_command() verbatim keeps the external run byte-identical to the in-VM run (so the ##VEP= header
# check passes on import). Every server path in the command is rewritten to a config.yaml placeholder so the
# compute box can install VEP + annotation data at different paths.
# --------------------------------------------------------------------------------------------------------

def _vep_command_config_roots() -> list[tuple[str, str]]:
    """ (config_key, server_path) roots used to rewrite absolute server paths in the VEP command into
        config.yaml placeholders. Everything VEP reads lives under these; annotation data + the fasta sit
        under annotation_base_dir. Longest server path first so nested dirs (cache/plugins/code) win over
        the base dir when both are prefixes of a token. """
    roots = [
        ("vep_code_dir", settings.ANNOTATION_VEP_CODE_DIR),
        ("vep_cache_dir", settings.ANNOTATION_VEP_CACHE_DIR),
        ("vep_plugins_dir", settings.ANNOTATION_VEP_PLUGINS_DIR),
        ("annotation_base_dir", settings.ANNOTATION_BASE_DIR),
    ]
    if settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT:
        roots.append(("perlbrew_runner", settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT))
    roots = [(key, value) for key, value in roots if value]
    roots.sort(key=lambda kv: len(kv[1]), reverse=True)
    return roots


def _templatize_scalar(template: list[str], flag: str, config_key: str) -> list[str]:
    """ Make a scalar VEP flag (e.g. --fork/--buffer_size) overridable from config: replace its value with a
        {config_key} placeholder, injecting the flag if get_vep_command() omitted it (fork==1). """
    template = list(template)
    placeholder = "{" + config_key + "}"
    if flag in template:
        index = template.index(flag)
        if index + 1 < len(template):
            template[index + 1] = placeholder
    else:
        template.extend([flag, placeholder])
    return template


def build_vep_command_template(variant_annotation_version: VariantAnnotationVersion,
                               pipeline_type=VariantAnnotationPipelineType.STANDARD) -> list[str]:
    """ The real VEP command (get_vep_command) with input/output as {input}/{output} and every server path
        rewritten to a {config-key} placeholder - rendered against config.yaml on the compute box (#1568). """
    genome_build = variant_annotation_version.genome_build
    annotation_consortium = variant_annotation_version.annotation_consortium
    cmd = get_vep_command("{input}", "{output}", genome_build, annotation_consortium, pipeline_type,
                          variant_annotation_version=variant_annotation_version)

    roots = _vep_command_config_roots()
    template = []
    for token in cmd:
        for config_key, server_path in roots:
            if server_path in token:
                token = token.replace(server_path, "{" + config_key + "}")
        template.append(token)

    template = _templatize_scalar(template, "--fork", "vep_fork")
    template = _templatize_scalar(template, "--buffer_size", "vep_buffer_size")
    return template


def build_snakemake_config(pipeline_type=VariantAnnotationPipelineType.STANDARD) -> dict:
    """ config.yaml defaults - server paths the operator edits to point at the compute box (#1568). """
    config = {
        "dump_dir": ".",
        "output_dir": "annotated",
        "vep_code_dir": settings.ANNOTATION_VEP_CODE_DIR,
        "vep_cache_dir": settings.ANNOTATION_VEP_CACHE_DIR,
        "vep_plugins_dir": settings.ANNOTATION_VEP_PLUGINS_DIR,
        "annotation_base_dir": settings.ANNOTATION_BASE_DIR,
        "vep_fork": settings.ANNOTATION_VEP_FORK or 1,
        "vep_buffer_size": settings.ANNOTATION_VEP_BUFFER_SIZE.get(pipeline_type) or 1000,
    }
    if settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT:
        config["perlbrew_runner"] = settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT
    return config


def _yaml_kv(key: str, value) -> str:
    if isinstance(value, int):
        return f"{key}: {value}"
    return f'{key}: "{value}"'


def _render_config_yaml(config: dict) -> str:
    lines = [
        "# VariantGrid external annotation config (#1568) - generated by `annotation_external --dump`.",
        "#",
        "# Edit the paths below to point at THIS compute box's VEP install + annotation data, then run:",
        "#     snakemake --cores <N>",
        "#",
        "# The directory tree under these roots must mirror the VariantGrid server's layout (VEP cache,",
        "# plugins, fasta and plugin data files all live under annotation_base_dir).",
        "",
        "# Dumped input VCFs + *.meta.json, and where to write annotated output (relative to this file or absolute):",
        _yaml_kv("dump_dir", config["dump_dir"]),
        _yaml_kv("output_dir", config["output_dir"]),
        "",
        "# VEP install + annotation data roots (server defaults shown - change to local paths):",
        _yaml_kv("vep_code_dir", config["vep_code_dir"]),
        _yaml_kv("vep_cache_dir", config["vep_cache_dir"]),
        _yaml_kv("vep_plugins_dir", config["vep_plugins_dir"]),
        _yaml_kv("annotation_base_dir", config["annotation_base_dir"]),
    ]
    if "perlbrew_runner" in config:
        lines.append(_yaml_kv("perlbrew_runner", config["perlbrew_runner"]))
    lines += [
        "",
        "# Compute knobs - tune for this box (more forks = more cores/RAM):",
        _yaml_kv("vep_fork", config["vep_fork"]),
        _yaml_kv("vep_buffer_size", config["vep_buffer_size"]),
        "",
    ]
    return "\n".join(lines)


# Self-contained Snakefile template (no VariantGrid imports so it runs standalone on the compute box). The
# VEP command template + annotated suffix are substituted for the __TEMPLATE_JSON__/__ANNOTATED_SUFFIX__
# sentinels at generation time; {input}/{output}/{config-key} braces are resolved by render_vep_command()
# at Snakemake runtime.
_SNAKEFILE_TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "templates_external", "Snakefile.template")


def render_snakefile(vep_command_template: list[str]) -> str:
    with open(_SNAKEFILE_TEMPLATE_PATH) as f:
        template = f.read()
    return (template
            .replace("__TEMPLATE_JSON__", json.dumps(vep_command_template, indent=4))
            .replace("__ANNOTATED_SUFFIX__", ANNOTATED_VCF_SUFFIX))


def write_snakemake_bundle(output_dir: str,
                           variant_annotation_version: VariantAnnotationVersion,
                           pipeline_type=VariantAnnotationPipelineType.STANDARD) -> tuple[str, str]:
    """ Write Snakefile + config.yaml into output_dir so the dumped VCFs can be annotated off-VM (#1568). """
    os.makedirs(output_dir, exist_ok=True)
    template = build_vep_command_template(variant_annotation_version, pipeline_type=pipeline_type)
    config = build_snakemake_config(pipeline_type=pipeline_type)

    config_path = os.path.join(output_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(_render_config_yaml(config))

    snakefile_path = os.path.join(output_dir, "Snakefile")
    with open(snakefile_path, "w") as f:
        f.write(render_snakefile(template))

    logging.info("Wrote Snakemake bundle: %s + %s", snakefile_path, config_path)
    return snakefile_path, config_path


# --------------------------------------------------------------------------------------------------------
# Import mode (#1568 §4.4): match annotated VCFs back to local EXTERNAL_DUMP_COMPLETED runs by version
# identity + range coordinate strings, run the §6a alignment pre-flight, then reuse the existing upload-only
# path. Per-run failure (bad match / §6a mismatch) never aborts the whole import - the run is skipped and
# falls back to a normal VEP run.
# --------------------------------------------------------------------------------------------------------

def _matchable_identity(identity: dict) -> dict:
    """ Version identity with the settings-snapshot fields excluded: distance/gencode_subset are re-derived
        from live settings and not in the ##VEP= header, so matching on them would produce false misses when
        a setting drifted after VAV creation (mirrors _vep_check_version_match). """
    return {key: value for key, value in identity.items() if key not in ("distance", "gencode_subset")}


def _expected_matchable_identity(meta: dict) -> dict:
    return _matchable_identity({key: value for key, value in meta["variant_annotation_version"].items()
                                if key != "pk"})


def _identity_diff(expected: dict, actual: dict) -> dict:
    """ field -> (dump_value, local_value) for every key that differs or is absent on one side. """
    diff = {}
    for key in set(expected) | set(actual):
        dump_value = expected.get(key, "<absent>")
        local_value = actual.get(key, "<absent>")
        if dump_value != local_value:
            diff[key] = (dump_value, local_value)
    return diff


def find_matching_variant_annotation_version(meta: dict) -> Optional[VariantAnnotationVersion]:
    """ Local VariantAnnotationVersion whose identity equals the annotated file's meta (#1568 §4.4). """
    genome_build = GenomeBuild.get_name_or_alias(meta["genome_build"])
    annotation_consortium = meta["annotation_consortium"]
    expected = _expected_matchable_identity(meta)
    for variant_annotation_version in VariantAnnotationVersion.objects.filter(
            genome_build=genome_build, annotation_consortium=annotation_consortium):
        if _matchable_identity(variant_annotation_version_identity(variant_annotation_version)) == expected:
            return variant_annotation_version
    return None


def explain_unmatched_variant_annotation_version(meta: dict) -> str:
    """ Human-readable reason no local VariantAnnotationVersion matched the file's version identity (#1568):
        the closest local VAV (fewest differing fields) plus the dump-vs-local values of each field that
        differs, so an operator can see exactly which annotation version is missing/misconfigured. """
    genome_build = GenomeBuild.get_name_or_alias(meta["genome_build"])
    annotation_consortium = meta["annotation_consortium"]
    expected = _expected_matchable_identity(meta)
    best = None  # (variant_annotation_version, diff)
    for variant_annotation_version in VariantAnnotationVersion.objects.filter(
            genome_build=genome_build, annotation_consortium=annotation_consortium):
        actual = _matchable_identity(variant_annotation_version_identity(variant_annotation_version))
        diff = _identity_diff(expected, actual)
        if best is None or len(diff) < len(best[1]):
            best = (variant_annotation_version, diff)
    if best is None:
        return f"no local VariantAnnotationVersion exists for {genome_build}/{annotation_consortium}"
    variant_annotation_version, diff = best
    fields = "; ".join(f"{key}: dump={dump_value!r} local={local_value!r}"
                       for key, (dump_value, local_value) in sorted(diff.items()))
    return (f"closest local VAV {variant_annotation_version.pk} ({variant_annotation_version}) "
            f"differs on {len(diff)} field(s): {fields}")


def find_matching_annotation_run(meta: dict, variant_annotation_version: VariantAnnotationVersion,
                                 pipeline_type=VariantAnnotationPipelineType.STANDARD) -> Optional[AnnotationRun]:
    """ Local external run still awaiting annotation whose range-lock min/max coordinate strings exactly equal
        the meta's (#1568 §3). Range locks are disjoint per version so at most one matches; annotation_run_pk
        is NOT used (it differs between a DB and its clone). """
    meta_range = meta["range"]
    candidates = AnnotationRun.objects.filter(
        annotation_range_lock__version=variant_annotation_version,
        pipeline_type=pipeline_type,
        external=True,
        vcf_annotated_filename__isnull=True,
    )
    for annotation_run in candidates:
        annotation_range_lock = annotation_run.annotation_range_lock
        if (str(annotation_range_lock.min_variant.coordinate) == meta_range["min_variant"]
                and str(annotation_range_lock.max_variant.coordinate) == meta_range["max_variant"]):
            return annotation_run
    return None


def explain_unmatched_annotation_run(meta: dict, variant_annotation_version: VariantAnnotationVersion,
                                     pipeline_type=VariantAnnotationPipelineType.STANDARD) -> str:
    """ Human-readable reason no external run awaiting annotation matched the file's range (#1568).
        Distinguishes the benign cases (the in-VM pipeline already annotated the range, or is still
        chewing on it, or it was already imported) from the real problem (no local range lock has this
        range at all - the file was produced against a different/diverged database), so an operator can
        tell a "nothing to do" skip apart from a genuine mismatch. """
    meta_range = meta["range"]
    min_coordinate, max_coordinate = meta_range["min_variant"], meta_range["max_variant"]
    matching_locks = [
        annotation_range_lock
        for annotation_range_lock in AnnotationRangeLock.objects.filter(version=variant_annotation_version)
        .select_related("min_variant", "max_variant")
        if str(annotation_range_lock.min_variant.coordinate) == min_coordinate
        and str(annotation_range_lock.max_variant.coordinate) == max_coordinate
    ]
    if not matching_locks:
        return (f"no local range lock for {variant_annotation_version} has range "
                f"{min_coordinate}..{max_coordinate} - range boundaries differ (produced against a "
                f"different/diverged database?)")

    runs = list(AnnotationRun.objects.filter(annotation_range_lock__in=matching_locks))

    def describe(annotation_run: AnnotationRun) -> str:
        state = "annotated" if annotation_run.vcf_annotated_filename else f"status={annotation_run.status}"
        return (f"run {annotation_run.pk} (external={annotation_run.external}, "
                f"pipeline_type={annotation_run.pipeline_type}, {state})")

    same_pipeline = [r for r in runs if r.pipeline_type == pipeline_type]
    for annotation_run in same_pipeline:
        if annotation_run.external and annotation_run.vcf_annotated_filename:
            return f"already imported - external {describe(annotation_run)}"
    for annotation_run in same_pipeline:
        if not annotation_run.external and annotation_run.vcf_annotated_filename:
            return (f"range already annotated by the in-VM pipeline - {describe(annotation_run)}; "
                    f"nothing to import")
    for annotation_run in same_pipeline:
        if not annotation_run.external:
            return f"range is on the in-VM pipeline, not offloaded - {describe(annotation_run)}"
    return ("no matching-pipeline external run awaiting annotation for this range; found "
            + ", ".join(describe(annotation_run) for annotation_run in runs))


def find_annotated_vcf(meta_path: str) -> Optional[str]:
    """ The annotated VCF sitting alongside a .meta.json sidecar (written by the Snakemake bundle) - None if
        it has not been annotated/copied back yet (#1568). """
    directory = os.path.dirname(meta_path)
    stem = os.path.basename(meta_path)[:-len(".meta.json")]
    for suffix in (ANNOTATED_VCF_SUFFIX, ".vep_annotated.vcf"):
        candidate = os.path.join(directory, stem + suffix)
        if os.path.exists(candidate):
            return candidate
    return None


def _import_annotated_annotation_run(annotation_run: AnnotationRun, annotated_vcf: str):
    """ Copy the annotated VCF into ANNOTATION_VCF_DUMP_DIR and hand the run to the normal single-authority
        dispatcher as a resume upload-only run (#1568). The off-VM VEP step is done, so the run rejoins the
        in-VM pipeline past VEP: we clear `external` (the dispatcher and its lease/reclaim system filter
        external=False, so an external run is invisible to them) and stamp the annotation start/end dates so
        get_status() -> ANNOTATION_COMPLETED - the dispatcher's priority-0 upload-only lane
        (_dispatchable_runs_qs resume clause). annotate_variants then skips dump+VEP straight to
        import_vcf_annotations, which runs the ##VEP= header check + inserts rows, setting upload_* -> FINISHED.

        We deliberately do NOT apply_async the upload ourselves: a raw queued celery job bypasses the
        dispatcher's capacity accounting + lease/reclaim (#2667/#1646), so it competes with the in-VM pipeline
        for worker slots and is lost for good on a worker restart. Routing through the dispatcher means a
        stalled upload is reclaimed and re-launched like any other run. """
    dest = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, os.path.basename(annotated_vcf))
    mk_path_for_file(dest)
    if os.path.abspath(annotated_vcf) != os.path.abspath(dest):
        shutil.copy(annotated_vcf, dest)

    now = timezone.now()
    annotation_run.vcf_annotated_filename = dest
    annotation_run.external = False  # VEP done off-VM - rejoin the in-VM pipeline for the DB upload
    annotation_run.annotation_start = now
    annotation_run.annotation_end = now
    annotation_run.upload_start = None
    annotation_run.upload_end = None
    annotation_run.error_exception = None
    annotation_run.task_id = None
    annotation_run.leased_by = None
    annotation_run.lease_expires = None
    annotation_run.save()  # get_status() -> ANNOTATION_COMPLETED (dispatcher upload-only lane)


def import_external_annotation_runs(genome_build: GenomeBuild, input_dir: str,
                                    pipeline_type=VariantAnnotationPipelineType.STANDARD,
                                    dry_run: bool = False, emit=None) -> dict:
    """ Match every annotated VCF in input_dir back to a local run and import it (#1568 §4.4). Returns a
        report of categorised outcomes; a per-file failure marks only that run and continues.

        `emit`, if given, is called emit(category, message) the moment each file's outcome is decided (the
        same category/message that lands in the report), so a caller can stream per-file progress + reasons
        live rather than waiting for the whole - potentially very long - run to finish. """
    report = {"imported": [], "matched": [], "unmatched": [], "missing_annotated": [], "id_mismatch": []}
    dispatched_vav_ids = set()

    def record(category: str, message: str):
        report[category].append(message)
        if emit is not None:
            emit(category, message)

    meta_paths = sorted(glob.glob(os.path.join(input_dir, "*.meta.json")))
    total = len(meta_paths)
    logging.info("External annotation import: %d meta file(s) in %s%s",
                 total, input_dir, " (dry-run)" if dry_run else "")
    last_progress = time.monotonic()
    for processed, meta_path in enumerate(meta_paths, start=1):
        now = time.monotonic()
        if now - last_progress >= IMPORT_PROGRESS_INTERVAL_SECONDS:
            logging.info(_import_progress_line(processed - 1, total, report))
            last_progress = now
        meta = parse_dump_metadata(meta_path)
        if meta.get("genome_build") != genome_build.name or meta.get("pipeline_type") != pipeline_type:
            continue

        variant_annotation_version = find_matching_variant_annotation_version(meta)
        if variant_annotation_version is None:
            record("unmatched",
                   f"{os.path.basename(meta_path)}: no local VariantAnnotationVersion matches version "
                   f"identity - {explain_unmatched_variant_annotation_version(meta)}")
            continue

        annotation_run = find_matching_annotation_run(meta, variant_annotation_version, pipeline_type)
        if annotation_run is None:
            record("unmatched",
                   f"{os.path.basename(meta_path)}: matched VAV {variant_annotation_version.pk} but "
                   f"{explain_unmatched_annotation_run(meta, variant_annotation_version, pipeline_type)}")
            continue

        annotated_vcf = find_annotated_vcf(meta_path)
        if annotated_vcf is None:
            record("missing_annotated",
                   f"{os.path.basename(meta_path)}: no annotated VCF alongside (not yet annotated?)")
            continue

        try:
            verify_annotated_vcf_variant_ids(annotation_run, meta)
        except VariantIdAlignmentError as e:
            annotation_run.error_exception = str(e)
            annotation_run.save()
            record("id_mismatch", str(e))
            continue

        if dry_run:
            record("matched",
                   f"{os.path.basename(annotated_vcf)} -> AnnotationRun {annotation_run.pk} "
                   f"({variant_annotation_version})")
            continue

        _import_annotated_annotation_run(annotation_run, annotated_vcf)
        dispatched_vav_ids.add(variant_annotation_version.pk)
        record("imported",
               f"{os.path.basename(annotated_vcf)} -> AnnotationRun {annotation_run.pk}")

    # Kick the single-authority dispatcher for each version we fed upload-only runs into. Offloaded runs live
    # on a NEW VariantAnnotationVersion, which the beat safety-net (ACTIVE-only) never sweeps - so without this
    # initial kick nothing would launch them. Once launched, each completion re-kicks the dispatcher for the
    # same version (annotate_variants._trigger_dispatch), draining the rest as worker slots free up.
    for vav_id in sorted(dispatched_vav_ids):
        dispatch_annotation_runs.si(vav_id).apply_async()

    logging.info(_import_progress_line(total, total, report))
    return report
