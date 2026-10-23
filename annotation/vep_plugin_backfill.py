"""
Switching on UTRAnnotator, AVI and SpliceVault (#2426) and the AlphaGenome splicing custom VCF (#2433) for
VariantAnnotationVersions already in use, without a new version: stamp their basename pins from live settings so the
VEP lane doesn't halt on VEPVersionMismatchError, flag the versions whose existing variants need a backfill, and write
the backfill steps.

Entry points: stamp_plugin_pins (the migration and `manage.py enable_vep_plugin_backfill`), backfill_steps /
backfill_note (the ManualOperation note and the command's output). Functions take the model classes, so the migration
passes its historical ones; nothing here imports annotation.models.
"""
from dataclasses import dataclass, field
from typing import Optional, Union

from django.conf import settings
from django.db.models import Max

from annotation import vep_columns
from annotation.enums import VEPCustom, VEPPlugin
from annotation.vep_config import vep_component_version_kwargs

VEPSource = Union[VEPPlugin, VEPCustom]

# VariantAnnotationVersion pin field -> the plugin or custom VCF whose data file it pins
SOURCE_PINS: dict[str, VEPSource] = {
    "utr_annotator": VEPPlugin.UTR_ANNOTATOR,
    "avi": VEPPlugin.AVI,
    "splicevault": VEPPlugin.SPLICEVAULT,
    "alphagenome_splicing": VEPCustom.ALPHAGENOME_SPLICING,
}
BACKFILL_FLAG = "backfilled_utr_avi_splicevault"
ACTIVE_STATUSES = ("NEW", "ACTIVE")  # HISTORICAL versions are never looked up by config again


@dataclass
class PinStamp:
    variant_annotation_version: object  # carries the stamped values, saved or (dry run) not
    changed_pins: dict[str, tuple] = field(default_factory=dict)  # field -> (old, new)

    @property
    def needs_backfill(self) -> bool:
        return not getattr(self.variant_annotation_version, BACKFILL_FLAG)


def source_pin_fields(variant_annotation_version_model) -> list[str]:
    """ The pins this model has - an earlier migration's historical model lacks the later sources' """
    model_fields = {f.name for f in variant_annotation_version_model._meta.get_fields()}
    return [pin_field for pin_field in SOURCE_PINS if pin_field in model_fields]


def source_columns(vav, vep_source: VEPSource) -> list:
    source_kwarg = "vep_custom" if isinstance(vep_source, VEPCustom) else "vep_plugin"
    return vep_columns.filter_for(**{source_kwarg: vep_source}, genome_build_name=vav.genome_build_id,
                                  columns_version=vav.columns_version, vep_version=vav.vep)


def stamp_plugin_pins(variant_annotation_version_model, annotation_range_lock_model,
                      genome_build_name: Optional[str] = None, dry_run: bool = False) -> list[PinStamp]:
    """ Sets each NEW/ACTIVE version's pins to what live settings pass VEP. A pin that changes to a file whose
        columns cover the version, on a version with annotated variants, sets the backfill flag False - those
        variants were annotated without it. An unchanged pin leaves the flag alone, so a re-run is a no-op. """
    qs = variant_annotation_version_model.objects.filter(status__in=ACTIVE_STATUSES)
    if genome_build_name:
        qs = qs.filter(genome_build_id=genome_build_name)

    pin_fields = source_pin_fields(variant_annotation_version_model)
    stamps = []
    for vav in qs.order_by("genome_build_id", "pk"):
        build_settings = settings.ANNOTATION.get(vav.genome_build_id)
        if build_settings is None:
            continue  # Build no longer configured on this deployment
        pins = vep_component_version_kwargs(build_settings)
        stamp = PinStamp(variant_annotation_version=vav)
        newly_covered = False
        for pin_field in pin_fields:
            old, new = getattr(vav, pin_field), pins[pin_field]
            if old != new:
                stamp.changed_pins[pin_field] = (old, new)
                setattr(vav, pin_field, new)
                newly_covered |= bool(new and source_columns(vav, SOURCE_PINS[pin_field]))
        if newly_covered and annotation_range_lock_model.objects.filter(version=vav).exists():
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


def backfill_steps(vav, annotation_range_lock_model) -> list[str]:
    """ One version's commands, bounded by the last variant it had locked for annotation when the plugins were
        switched on - anything after that is annotated with them """
    build = vav.genome_build_id
    columns_by_source = backfill_source_columns(vav)
    columns = {column: csq_field for c in columns_by_source.values() for column, csq_field in c.items()}
    max_variant_id = annotation_range_lock_model.objects.filter(version=vav).aggregate(Max("max_variant_id"))
    common = f"--genome-build {build} --variant-annotation-version {vav.pk} --columns {','.join(columns)}"
    csq_args = " ".join(f"--csq {column}={csq_field}" for column, csq_field in columns.items())
    csq_fields = ",".join(["Feature", "PICK", *columns.values()])
    sources = ", ".join(vep_source.label for vep_source in columns_by_source)
    return [
        f"{build} (VariantAnnotationVersion {vav.pk}):",
        f"1. python3 manage.py annotation_backfill_columns --dump {common} "
        f"--max-variant-id {max_variant_id['max_variant_id__max']} --output vep_plugins_{build}.vcf.gz",
        f"2. Annotate it with VEP: this build's STANDARD 'Show VEP command' on /annotation/annotation_versions, with "
        f"-i vep_plugins_{build}.vcf.gz -o vep_plugins_{build}.vep_annotated.vcf.gz, every --plugin and --custom "
        f"except {sources} removed, and --fields {csq_fields}. Keep the other flags, so the PICK'd transcript matches",
        f"3. python3 manage.py annotation_backfill_columns --import {common} {csq_args} "
        f"--vcf vep_plugins_{build}.vep_annotated.vcf.gz --dry-run",
        f"4. The same without --dry-run, plus --mark-backfilled {BACKFILL_FLAG}",
    ]


def backfill_note(variant_annotation_version_model, annotation_range_lock_model,
                  genome_build_name: Optional[str] = None) -> str:
    """ Steps for every version still flagged for the backfill """
    qs = variant_annotation_version_model.objects.filter(**{BACKFILL_FLAG: False})
    if genome_build_name:
        qs = qs.filter(genome_build_id=genome_build_name)
    lines = ["Fill in UTRAnnotator, AVI, SpliceVault and AlphaGenome splicing annotation of variants annotated "
             "before they were switched on (#2426, #2433), with one VEP run per version. Full recipe in claude/research/annotation.md, 'A new column "
             "on an existing VariantAnnotationVersion'. For each version:"]
    for vav in qs.order_by("genome_build_id", "pk"):
        lines.append("")
        lines.extend(backfill_steps(vav, annotation_range_lock_model))
    return "\n".join(lines)
