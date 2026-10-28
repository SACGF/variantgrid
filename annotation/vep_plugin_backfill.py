"""
Switching on UTRAnnotator, AVI and SpliceVault (#2426) and the AlphaGenome splicing custom VCF (#2433) for
VariantAnnotationVersions already in use, without a new version: stamp their basename pins from live settings so the
VEP lane doesn't halt on VEPVersionMismatchError, flag the versions whose existing variants need a backfill, and write
the backfill steps.

Entry points: stamp_plugin_pins (the migration and `manage.py enable_vep_plugin_backfill`), backfill_steps
(the command's output, with the backfill's VEP command line built by backfill_vep_command), backfill_script (the same
steps as a bash script, the command's `--bash-script`), backfill_note (the ManualOperation note, pointing at the script),
enable_plugins_fix / vep_mismatch_fix (what the mismatch warnings tell the admin to run). Functions take the model
classes, so the migration passes its historical ones; nothing here imports annotation.models.
"""
import os
import shlex
from dataclasses import dataclass, field
from typing import Optional, Union

from django.conf import settings
from django.db.models import Max, Q

from annotation import vep_columns
from annotation.enums import VariantAnnotationPipelineType, VEPCustom, VEPPlugin
from annotation.vep_annotation import get_vep_command, get_vep_skipped_variants_filename
from annotation.vep_config import vep_component_version_kwargs
from snpdb.models.models_genome import GenomeBuild

VEPSource = Union[VEPPlugin, VEPCustom]

# VariantAnnotationVersion pin field -> the plugin or custom VCF whose data file it pins
SOURCE_PINS: dict[str, VEPSource] = {
    "utr_annotator": VEPPlugin.UTR_ANNOTATOR,
    "avi": VEPPlugin.AVI,
    "splicevault": VEPPlugin.SPLICEVAULT,
    "alphagenome_splicing": VEPCustom.ALPHAGENOME_SPLICING,
}
BACKFILL_FLAG = "backfilled_utr_avi_splicevault"
# Variants annotated by Ensembl's UTRAnnotator, which dropped all but one of a variant's effects at random (#2499): their
# utr_annotator columns are re-imported from a UTRAnnotator-only pass. A version still waiting for the plugin backfill
# above gets both in that one pass, over every variant rather than only those annotated before the switch-on
UTR_EFFECTS_FLAG = "backfilled_utr_effects"
UTR_EFFECTS_COLUMNS = ("utr_annotator_consequence", "utr_annotator_annotation")
# The utr_annotator_* and splicevault_* columns are transcript-level - the variant page's transcript switcher and
# evidence key autopopulate read them off the transcript rows - but the plugin backfill imported the PICK'd entry into
# VariantAnnotation only, leaving the transcript rows null. They are filled from a UTRAnnotator + SpliceVault pass over
# every variant, each CSQ entry imported to the row of its transcript
TRANSCRIPTS_FLAG = "backfilled_utr_splicevault_transcripts"
TRANSCRIPT_SOURCES = (VEPPlugin.UTR_ANNOTATOR, VEPPlugin.SPLICEVAULT)
ENABLE_COMMAND = "python3 manage.py enable_vep_plugin_backfill"
DOWNLOAD_COMMAND = "python3 manage.py download_annotation_data"
ACTIVE_STATUSES = ("NEW", "ACTIVE")  # HISTORICAL versions are never looked up by config again
# Flags that only add CSQ fields the backfill doesn't import. Dropping them leaves the PICK'd entry the same and runs
# VEP ~1.5x faster
_OUTPUT_ONLY_FLAGS = {"--hgvs", "--check_existing", "--af", "--pubmed", "--domains", "--protein", "--numbers"}
_OUTPUT_ONLY_FLAGS_WITH_VALUE = {"--sift"}
# The pipeline's buffer suits its small batches; one run of every variant goes faster with each fork given more
_BACKFILL_BUFFER_SIZE = 5000


@dataclass
class PinStamp:
    variant_annotation_version: object  # carries the stamped values, saved or (dry run) not
    changed_pins: dict[str, tuple] = field(default_factory=dict)  # field -> (old, new)

    @property
    def needs_backfill(self) -> bool:
        return bool(pending_backfill_flags(self.variant_annotation_version))


@dataclass
class BackfillStep:
    command: str  # one shell command, run from the VariantGrid checkout
    comment: Optional[str] = None  # why it's there, where the command alone doesn't say


def source_pin_fields(variant_annotation_version_model) -> list[str]:
    """ The pins this model has - an earlier migration's historical model lacks the later sources' """
    model_fields = {f.name for f in variant_annotation_version_model._meta.get_fields()}
    return [pin_field for pin_field in SOURCE_PINS if pin_field in model_fields]


def source_columns(vav, vep_source: VEPSource) -> list:
    source_kwarg = "vep_custom" if isinstance(vep_source, VEPCustom) else "vep_plugin"
    return vep_columns.filter_for(**{source_kwarg: vep_source}, genome_build_name=vav.genome_build_id,
                                  columns_version=vav.columns_version, vep_version=vav.vep)


def settings_pin_changes(vav) -> dict[str, tuple]:
    """ pin field -> (version's pin, live settings' pin) for each source whose settings file the version doesn't pin.
        Reads settings only, so a web node gives the same answer as an annotation node """
    build_settings = settings.ANNOTATION.get(vav.genome_build_id)
    if build_settings is None:
        return {}
    pins = vep_component_version_kwargs(build_settings)
    return {pin_field: (getattr(vav, pin_field), pins[pin_field]) for pin_field in source_pin_fields(type(vav))
            if getattr(vav, pin_field) != pins[pin_field]}


def _enable_fix(pin_fields) -> str:
    return (f"{', '.join(pin_fields)} set in settings but not pinned on this version, so every VEP run against it "
            f"fails the version check. Stop annotation_workers, run '{ENABLE_COMMAND} --dry-run' then "
            f"'{ENABLE_COMMAND}' to stamp the pins, restart the celery workers and gunicorn, and run the backfill "
            f"steps it prints for variants already annotated (#2426, #2433)")


def enable_plugins_fix(vav) -> Optional[str]:
    """ What to run when settings switch on a plugin / custom VCF this version doesn't pin, else None """
    if vav.status in ACTIVE_STATUSES and (changes := settings_pin_changes(vav)):
        return _enable_fix(changes)
    return None


def vep_mismatch_fix(vav, vep_vav_kwargs: dict) -> str:
    """ What to run when vav differs from vep_vav_kwargs (get_vep_variant_annotation_version_kwargs): stamping
        the plugin pins when only those differ, otherwise a new version """
    mismatched = [k for k, v in vep_vav_kwargs.items() if getattr(vav, k) != v]
    if vav.status in ACTIVE_STATUSES and mismatched and set(mismatched) <= SOURCE_PINS.keys():
        return _enable_fix(mismatched)
    return ("Upgrade via 'python3 manage.py create_new_variant_annotation_version "
            f"--genome-build={vav.genome_build_id}'")


def stamp_plugin_pins(variant_annotation_version_model, annotation_range_lock_model,
                      genome_build_name: Optional[str] = None, dry_run: bool = False) -> list[PinStamp]:
    """ Sets each NEW/ACTIVE version's pins to what live settings pass VEP, so none fails the version check. A pin
        that changes to a file whose columns cover the version, on the build's newest version with annotated
        variants, sets the backfill flag False - those variants were annotated without it. Older versions are on
        their way out, so aren't worth a backfill. An unchanged pin leaves the flag alone, so a re-run is a no-op. """
    qs = variant_annotation_version_model.objects.filter(status__in=ACTIVE_STATUSES)
    if genome_build_name:
        qs = qs.filter(genome_build_id=genome_build_name)
    vavs = list(qs.order_by("genome_build_id", "annotation_date", "pk"))
    newest_pk_by_build = {vav.genome_build_id: vav.pk for vav in vavs}

    stamps = []
    for vav in vavs:
        if vav.genome_build_id not in settings.ANNOTATION:
            continue  # Build no longer configured on this deployment
        stamp = PinStamp(variant_annotation_version=vav, changed_pins=settings_pin_changes(vav))
        newly_covered = False
        for pin_field, (_old, new) in stamp.changed_pins.items():
            setattr(vav, pin_field, new)
            newly_covered |= bool(new and source_columns(vav, SOURCE_PINS[pin_field]))
        is_newest = vav.pk == newest_pk_by_build[vav.genome_build_id]
        if newly_covered and is_newest and annotation_range_lock_model.objects.filter(version=vav).exists():
            setattr(vav, BACKFILL_FLAG, False)
        if stamp.changed_pins and not dry_run:
            vav.save(update_fields=[*stamp.changed_pins, BACKFILL_FLAG])
        stamps.append(stamp)
    return stamps


def backfill_source_columns(vav) -> dict[VEPSource, dict[str, str]]:
    """ VariantGrid column -> CSQ field, for each source this version pins whose columns it covers """
    columns_by_source = {}
    for pin_field in source_pin_fields(type(vav)):
        if getattr(vav, pin_field):
            vep_source = SOURCE_PINS[pin_field]
            if column_defs := source_columns(vav, vep_source):
                columns_by_source[vep_source] = {c.variant_grid_columns[0]: c.vep_info_field for c in column_defs}
    return columns_by_source


def backfill_fork() -> int:
    """ The pipeline usually runs several VEPs at once with ANNOTATION_VEP_FORK 1 - the backfill is one big run """
    return max(settings.ANNOTATION_VEP_FORK or 1, (os.cpu_count() or 2) - 1)


def backfill_vep_command(vav, vep_sources, csq_fields: list[str], input_vcf: str, output_vcf: str) -> list[str]:
    """ The pipeline's STANDARD command for this version, so the cache, pick order and PICK'd transcript match, with
        only the backfilled sources' plugins / custom files, their CSQ fields, and the output-only flags dropped """
    genome_build = GenomeBuild.get_name_or_alias(vav.genome_build_id)
    pipeline_cmd = get_vep_command(input_vcf, output_vcf, genome_build, vav.annotation_consortium,
                                   VariantAnnotationPipelineType.STANDARD, variant_annotation_version=vav)
    plugin_names = {s.label for s in vep_sources if isinstance(s, VEPPlugin)}
    custom_short_names = {f"short_name={s.label}" for s in vep_sources if isinstance(s, VEPCustom)}

    cmd = []
    args = iter(pipeline_cmd)
    for arg in args:
        if arg in _OUTPUT_ONLY_FLAGS:
            continue
        if arg in _OUTPUT_ONLY_FLAGS_WITH_VALUE or arg in ("--fork", "--buffer_size"):
            next(args)
            continue
        if arg == "--plugin":
            plugin_arg = next(args)
            if plugin_arg.split(",")[0] in plugin_names:
                cmd.extend([arg, plugin_arg])
            continue
        if arg == "--custom":
            custom_arg = next(args)
            if custom_short_names.intersection(custom_arg.split(",")):
                cmd.extend([arg, custom_arg])
            continue
        cmd.append(arg)
    cmd.extend(["--fork", str(backfill_fork()), "--buffer_size", str(_BACKFILL_BUFFER_SIZE),
                "--fields", ",".join(csq_fields)])
    return cmd


def pending_backfill_flags(vav) -> list[str]:
    """ The backfilled_* flags this version still has False, of the ones one VEP pass here can set """
    return [flag for flag in (BACKFILL_FLAG, UTR_EFFECTS_FLAG, TRANSCRIPTS_FLAG) if not getattr(vav, flag, True)]


def backfill_commands(vav, annotation_range_lock_model) -> list[BackfillStep]:
    """ One version's commands: one VEP pass and one import, whichever flags it has pending.

        The plugin backfill alone is bounded by the last variant locked when the plugins were switched on - anything
        after that is annotated with them. The transcript rows pass runs UTRAnnotator and SpliceVault over every
        variant (the rows of variants annotated before the switch-on are null whatever their representative row
        holds). The UTR effects pass alone runs UTRAnnotator over the variants that have a consequence, the only ones
        the old plugin could have mis-reported. A version needing the plugin backfill as well gets one unbounded pass
        with every source: re-importing the other sources' values for already annotated variants changes nothing. """
    flags = pending_backfill_flags(vav)
    if not flags:
        return []
    build = vav.genome_build_id
    columns_by_source = backfill_source_columns(vav)
    dump_args = ["--pipeline-type STANDARD"]
    if BACKFILL_FLAG in flags:
        if flags == [BACKFILL_FLAG]:
            max_variant_id = annotation_range_lock_model.objects.filter(version=vav).aggregate(Max("max_variant_id"))
            dump_args.append(f"--max-variant-id {max_variant_id['max_variant_id__max']}")
    elif TRANSCRIPTS_FLAG in flags:
        columns_by_source = {source: columns for source, columns in columns_by_source.items()
                             if source in TRANSCRIPT_SOURCES}
    else:
        columns_by_source = {source: {column: csq_field for column, csq_field in columns.items()
                                      if column in UTR_EFFECTS_COLUMNS}
                             for source, columns in columns_by_source.items() if source == VEPPlugin.UTR_ANNOTATOR}
        dump_args.append(f"--not-null {UTR_EFFECTS_COLUMNS[0]}")
    columns = {column: csq_field for c in columns_by_source.values() for column, csq_field in c.items()}
    common = f"--genome-build {build} --variant-annotation-version {vav.pk} --columns {','.join(columns)}"
    csq_args = " ".join(f"--csq {column}={csq_field}" for column, csq_field in columns.items())
    dump_vcf = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, f"vep_plugins_{build}.vcf.gz")
    annotated_vcf = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, f"vep_plugins_{build}.vep_annotated.vcf.gz")
    vep_cmd = backfill_vep_command(vav, columns_by_source, ["Feature", "PICK", *columns.values()],
                                   dump_vcf, annotated_vcf)
    import_cmd = (f"python3 manage.py annotation_backfill_columns --import {common} {csq_args} "
                  f"--vcf {annotated_vcf} --dump-vcf {dump_vcf}")
    steps = []
    if UTR_EFFECTS_FLAG in flags:
        # The pass is only a fix with our UTRAnnotator.pm installed over Ensembl's; the import checks before marking
        steps.append(BackfillStep(DOWNLOAD_COMMAND,
                                  "installs VariantGrid's UTRAnnotator.pm over the plugins dir's copy"))
    steps.extend([
        BackfillStep(f"python3 manage.py annotation_backfill_columns --dump {common} {' '.join(dump_args)} "
                     f"--output {dump_vcf}"),
        BackfillStep(shlex.join(vep_cmd)),
        BackfillStep(f"{import_cmd} --dry-run"),
        BackfillStep(f"{import_cmd} --mark-backfilled {','.join(flags)}"),
        BackfillStep(f"rm -f {dump_vcf} {annotated_vcf} {annotated_vcf}_warnings.txt "
                     f"{get_vep_skipped_variants_filename(annotated_vcf)}"),
    ])
    return steps


def _backfill_heading(vav) -> str:
    return f"{vav.genome_build_id} (VariantAnnotationVersion {vav.pk})"


def backfill_steps(vav, annotation_range_lock_model) -> list[str]:
    """ backfill_commands as numbered text under a heading, for the command's output """
    steps = backfill_commands(vav, annotation_range_lock_model)
    if not steps:
        return []
    lines = [f"{_backfill_heading(vav)}:"]
    for i, step in enumerate(steps, 1):
        comment = f"  ({step.comment})" if step.comment else ""
        lines.append(f"{i}. {step.command}{comment}")
    return lines


def backfill_script(vavs, annotation_range_lock_model) -> str:
    """ The backfill for these versions as one bash script: each version's steps in order, stopping at the first
        command that fails. The import --dry-run stays in: its 'would change N rows' is in the log, and it refuses
        VEP output that doesn't cover every dumped record before the real import gets to """
    lines = [
        "#!/bin/bash",
        f"# Generated by '{ENABLE_COMMAND} --bash-script'",
        "# Fills in the VEP plugin / custom VCF columns of variants annotated before they were switched on",
        "# (#2426, #2433), their UTRAnnotator / SpliceVault transcript rows, and re-annotates the 5' UTR effects",
        "# Ensembl's UTRAnnotator reported one of at random (#2499): one VEP pass per VariantAnnotationVersion.",
        "# Hours per version - run under screen or nohup.",
        "# Stops at the first failing command; a version whose import marked it backfilled is left out when the",
        "# script is generated again.",
        "set -euo pipefail",
        f"cd {shlex.quote(settings.BASE_DIR)}",
    ]
    versions = 0
    for vav in vavs:
        steps = backfill_commands(vav, annotation_range_lock_model)
        if not steps:
            continue
        versions += 1
        lines.extend(["", f"echo \"=== {_backfill_heading(vav)} $(date) ===\""])
        for step in steps:
            if step.comment:
                lines.append(f"# {step.comment}")
            lines.append(step.command)
    if not versions:
        lines.extend(["", "# No version is flagged for a backfill."])
    lines.extend(["", "echo \"=== Backfill complete $(date) ===\"", ""])
    return "\n".join(lines)


def versions_pending_backfill(variant_annotation_version_model, genome_build_name: Optional[str] = None):
    """ Every version still flagged for a backfill, in the order the steps run """
    pending = Q()
    for flag in (BACKFILL_FLAG, UTR_EFFECTS_FLAG, TRANSCRIPTS_FLAG):
        pending |= Q(**{flag: False})
    qs = variant_annotation_version_model.objects.filter(pending)
    if genome_build_name:
        qs = qs.filter(genome_build_id=genome_build_name)
    return qs.order_by("genome_build_id", "pk")


def backfill_note(variant_annotation_version_model, genome_build_name: Optional[str] = None) -> str:
    """ What the backfill is and which versions need it - the steps themselves come from the command, whose
        VEP command lines would bury the upgrader in output """
    versions = ", ".join(_backfill_heading(vav)
                         for vav in versions_pending_backfill(variant_annotation_version_model, genome_build_name))
    return ("Fill in UTRAnnotator, AVI, SpliceVault and AlphaGenome splicing annotation of variants annotated "
            "before they were switched on (#2426, #2433), the UTRAnnotator and SpliceVault transcript rows the "
            "earlier backfills left null, and re-annotate the 5' UTR effects Ensembl's UTRAnnotator "
            f"reported one of at random (#2499), with one VEP run per version: {versions}. "
            f"Write the steps as one script with '{ENABLE_COMMAND} --bash-script backfill.sh' and run it under "
            "screen or nohup (hours per version). Admin steps: "
            "https://github.com/SACGF/variantgrid/issues/2433#issuecomment-6393595562")
