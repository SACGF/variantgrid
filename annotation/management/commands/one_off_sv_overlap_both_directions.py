"""
Re-applies the gnomAD SV overlap rule to SV annotations made before SVOverlapProcessor checked both directions
(#2589). VEP's overlap_cutoff only requires the gnomAD SV to be covered, so a large SV took the lowest AF among
the small gnomAD SVs inside it. Every overlapping record is still on the row (gnomad_sv_overlap_name / _af /
_percent / _coords), so each VariantAnnotationVersion's SV rows are re-filtered here without VEP: a record below
sv_overlap_min_fraction in either direction leaves the multi-value fields, and when it was the record the gnomad_*
columns came from (the lowest AF), they are refilled from the lowest-AF survivor, read from the gnomAD SV VCF the
build is configured with, or blanked when none survive. A survivor missing from that VCF (another gnomAD release)
leaves its row as it was and is reported. One transaction per version; --dry-run reports what would change.
Registered as a deploy-time ManualOperation (annotation/migrations/0013_sv_overlap_both_directions.py).
"""
import logging
import os
from dataclasses import dataclass, field

import cyvcf2
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db import transaction

from annotation import vep_columns
from annotation.enums import VariantAnnotationPipelineType, VEPCustom
from annotation.models.models import VariantAnnotation
from annotation.models.models_annotation_run import AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_config import VEPConfig
from annotation.vep_field_formatters import VEP_SEPARATOR
from annotation.vep_plugin_display import (
    GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS,
    vep_multi_fields_to_list_of_dicts,
)
from library.genomics import Range, overlap_fraction, parse_gnomad_coord

STRUCTURAL_VARIANT = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


@dataclass
class RowRepair:
    variant_annotation: VariantAnnotation
    survivors: list[dict]
    refill_from: str | None = None  # gnomAD SV name whose INFO refills the gnomad_* columns


@dataclass
class VersionRepair:
    checked: int = 0
    pruned: int = 0  # a dropped record wasn't the one the gnomad_* columns came from
    refilled: int = 0  # it was: gnomad_* refilled from the lowest-AF survivor
    blanked: int = 0  # nothing survived
    missing: set[str] = field(default_factory=set)  # survivors not in the configured gnomAD SV VCF
    skipped: str | None = None

    def __str__(self):
        if self.skipped:
            return f"skipped: {self.skipped}"
        s = f"{self.checked} SV rows checked, {self.pruned} pruned, {self.refilled} refilled, {self.blanked} blanked"
        if self.missing:
            s += f", {len(self.missing)} left as they were (not in the gnomAD SV VCF: {sorted(self.missing)[:5]})"
        return s


def _gnomad_sv_defs(vav: VariantAnnotationVersion, vep_custom: VEPCustom):
    return vep_columns.filter_for(genome_build_name=vav.genome_build.name, pipeline_type=STRUCTURAL_VARIANT,
                                  columns_version=vav.columns_version, vep_version=vav.vep, vep_custom=vep_custom)


def _our_percent(ours: Range, gnomad: Range) -> int:
    return int(overlap_fraction(ours, gnomad) * 100)


def _passes_both_ways(ours: Range, record: dict, min_fraction: float) -> bool:
    _chrom, start, end = parse_gnomad_coord(record["gnomad_sv_overlap_coords"])
    gnomad = Range(start, end)
    return overlap_fraction(ours, gnomad) >= min_fraction and overlap_fraction(gnomad, ours) >= min_fraction


def find_row_repairs(vav: VariantAnnotationVersion, min_fraction: float) -> tuple[int, list[RowRepair]]:
    """ (SV rows checked, the rows with a record that fails the rule) """
    sv_run_ids = list(AnnotationRun.objects.filter(annotation_range_lock__version=vav,
                                                   pipeline_type=STRUCTURAL_VARIANT).values_list("pk", flat=True))
    checked = 0
    repairs = []
    if not sv_run_ids:
        return checked, repairs
    qs = VariantAnnotation.objects.filter(version=vav, annotation_run_id__in=sv_run_ids,
                                          gnomad_sv_overlap_coords__isnull=False)
    qs = qs.exclude(gnomad_sv_overlap_coords="").select_related("variant__locus")
    for va in qs.iterator():
        checked += 1
        data = {f: getattr(va, f) for f in GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS}
        records = vep_multi_fields_to_list_of_dicts(data, GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS)
        ours = Range(va.variant.locus.position, va.variant.end)
        survivors = [r for r in records if _passes_both_ways(ours, r, min_fraction)]
        if len(survivors) == len(records):
            continue
        for r in survivors:
            _chrom, start, end = parse_gnomad_coord(r["gnomad_sv_overlap_coords"])
            r["gnomad_sv_overlap_percent"] = str(_our_percent(ours, Range(start, end)))
        repair = RowRepair(va, survivors)
        # lowest_af is the only pick implemented, so the chosen record is the lowest AF of all of them
        chosen = min(records, key=lambda r: float(r["gnomad_sv_overlap_af"]))
        if survivors and chosen not in survivors:
            repair.refill_from = min(survivors, key=lambda r: float(r["gnomad_sv_overlap_af"]))["gnomad_sv_overlap_name"]
        repairs.append(repair)
    return checked, repairs


def read_gnomad_sv_records(vcf_path: str, names: set[str], info_fields: set[str]) -> dict[str, dict]:
    """ name -> INFO values, for the records in names (one pass, the file has no ID index) """
    records = {}
    for v in cyvcf2.VCF(vcf_path):
        if v.ID in names:
            records[v.ID] = {f: v.INFO.get(f) for f in info_fields}
            if len(records) == len(names):
                break
    return records


def repair_version(vav: VariantAnnotationVersion, dry_run: bool) -> VersionRepair:
    result = VersionRepair()
    if vav.sv_overlap_single_value_method != "lowest_af":
        result.skipped = f"sv_overlap_single_value_method {vav.sv_overlap_single_value_method!r} isn't lowest_af"
        return result
    min_fraction = vav.sv_overlap_min_fraction
    if min_fraction is None:
        min_fraction = settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION
    result.checked, repairs = find_row_repairs(vav, min_fraction)
    if not repairs:
        return result

    # gnomad_af etc. come from the chosen record's INFO fields; the multi-value fields from the survivors
    single_columns = {vgc: d.source_field
                      for d in _gnomad_sv_defs(vav, VEPCustom.GNOMAD_SV) if d.source_field
                      for vgc in d.variant_grid_columns if vgc not in GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS}
    all_columns = {vgc for custom in (VEPCustom.GNOMAD_SV, VEPCustom.GNOMAD_SV_NAME)
                   for d in _gnomad_sv_defs(vav, custom) for vgc in d.variant_grid_columns}

    gnomad_records = {}
    if refill_names := {r.refill_from for r in repairs if r.refill_from}:
        try:
            vcf_path = VEPConfig(vav.genome_build)["gnomad_sv"]
        except KeyError:
            vcf_path = None
        if vcf_path and os.path.exists(vcf_path):
            gnomad_records = read_gnomad_sv_records(vcf_path, refill_names, set(single_columns.values()))
        else:
            # A build that isn't enabled here has no annotation data downloaded
            logging.warning("%s: no gnomad_sv file for %s (%s), rows needing a refill are left as they were",
                            vav, vav.genome_build, vcf_path or "not configured")

    with transaction.atomic():
        for repair in repairs:
            if repair.refill_from and repair.refill_from not in gnomad_records:
                result.missing.add(repair.refill_from)
                continue
            if not repair.survivors:
                update = dict.fromkeys(all_columns)
                result.blanked += 1
            else:
                update = {f: VEP_SEPARATOR.join(r[f] for r in repair.survivors)
                          for f in GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS}
                if repair.refill_from:
                    info = gnomad_records[repair.refill_from]
                    update.update({vgc: info[source_field] for vgc, source_field in single_columns.items()})
                    result.refilled += 1
                else:
                    result.pruned += 1
            if not dry_run:
                VariantAnnotation.objects.filter(version=vav, pk=repair.variant_annotation.pk).update(**update)
    return result


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true", help="Report what would change without writing")
        parser.add_argument("--variant-annotation-version", type=int,
                            help="Only this version (default: every version whose data isn't archived)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        qs = VariantAnnotationVersion.objects.filter(data_archived_date__isnull=True).select_related("genome_build")
        if vav_id := options["variant_annotation_version"]:
            qs = qs.filter(pk=vav_id)
        for vav in qs.order_by("pk"):
            result = repair_version(vav, dry_run)
            prefix = "DRY RUN " if dry_run else ""
            self.stdout.write(f"{prefix}{vav}: {result}")
