"""
Re-pick the gnomAD-SV overlap record whose values were copied onto the gnomAD columns.

`SVOverlapProcessor._pick_record` compared VEP's raw string values when applying the 'lowest_af'
method, so wherever gnomAD-SV mixes scientific ('4.6e-05') and decimal ('0.006085') notation it
picked lexicographically and landed on the *most* common overlapping record. The chosen record's
values are what get written to gnomad_af / gnomad_ac / gnomad_an / gnomad_hom_alt /
gnomad_popmax_af and the per-population AFs - which is what the analysis PopulationNode filters
on - so affected SVs were being filtered out on an allele frequency belonging to the wrong record.

Only the '&'-joined gnomad_sv_overlap_* text fields survive per-record, so the correct record's
remaining columns are read back out of the configured gnomAD-SV VCF, matched on the stored name
and confirmed by its AF - which is what tells us it is the release the version was annotated with.

Entry point: manage.py fix_gnomad_sv_overlap_chosen_record [--dry-run] [--variant-annotation-version ID]
"""
import logging
from collections import Counter
from typing import Optional

import cyvcf2
from django.core.management.base import BaseCommand
from django.db import transaction

from annotation import vep_columns
from annotation.models import VariantAnnotation, VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType, VEPCustom
from annotation.vep_config import VEPConfig
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from library.genomics import parse_gnomad_coord
from snpdb.models import GenomeBuild

PROGRESS_INTERVAL = 1000


def _sv_column_targets(genome_build: GenomeBuild) -> dict[str, list[str]]:
    """ gnomAD-SV VCF INFO field -> VariantAnnotation columns holding the chosen record's value.
        The '&'-joined multi-value fields keep every overlapping record so are left alone. """
    targets = {}
    for column_def in vep_columns.filter_for(vep_custom=VEPCustom.GNOMAD_SV,
                                             genome_build_name=genome_build.name,
                                             pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
        if not column_def.source_field:
            continue  # VEP's own overlap output (coords/percent), not a source VCF INFO field
        columns = [vgc for vgc in column_def.variant_grid_columns
                   if vgc not in VariantAnnotation.GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS]
        if columns:
            targets[column_def.source_field] = columns
    return targets


def _raw_info(record: cyvcf2.Variant) -> dict[str, str]:
    """ cyvcf2 parses Float INFO through float32, widening '4.6e-05' to 4.600000102072954e-05.
        VEP passes the text straight through, so read the raw INFO column to store what
        re-annotating would have stored. """
    info_column = str(record).split("\t")[7]
    values = {}
    for entry in info_column.split(";"):
        key, _, value = entry.partition("=")
        values[key] = value
    return values


class GnomADSVRecordReader:
    """ Reads whole gnomAD-SV records back out of the VCF configured for a build, by name.
        The stored coords came from that same file, so they give us the tabix region. """

    def __init__(self, genome_build: GenomeBuild):
        self.vcf_filename = VEPConfig(genome_build)["gnomad_sv"]
        self.reader = cyvcf2.VCF(self.vcf_filename)
        self._cache = {}

    def get_info(self, name: str, coord: str) -> Optional[dict[str, str]]:
        cached = self._cache.get(name)
        if cached is None:
            chrom, start, _end = parse_gnomad_coord(coord)
            for record in self.reader(f"{chrom}:{start}-{start}"):
                if record.ID == name:
                    cached = _raw_info(record)
                    self._cache[name] = cached
                    break
        return cached


def _partition_qs(vav: VariantAnnotationVersion):
    """ Reads and writes go against this version's partition - the base table holds no rows """
    qs = get_queryset_with_transformer_hook(klass=VariantAnnotation)
    qs.add_sql_transformer(vav.sql_partition_transformer)
    return qs


def fix_variant_annotation_version(vav: VariantAnnotationVersion, dry_run: bool = False) -> Counter:
    results = Counter()

    targets = _sv_column_targets(vav.genome_build)
    try:
        reader = GnomADSVRecordReader(vav.genome_build)
    except KeyError:
        logging.info("%s has no gnomAD-SV configured for %s - skipping", vav, vav.genome_build)
        results["no_gnomad_sv_configured"] += 1
        return results

    partition_qs = _partition_qs(vav)
    qs = partition_qs.filter(gnomad_sv_overlap_af__contains="&")
    fields = ["variant_id", "gnomad_af", *VariantAnnotation.GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS]
    for va in qs.values(*fields).iterator():
        results["checked"] += 1

        af_values = [float(af) for af in va["gnomad_sv_overlap_af"].split("&")]
        lowest_af = min(af_values)
        if va["gnomad_af"] == lowest_af:
            results["already_correct"] += 1
            continue

        chosen_index = af_values.index(lowest_af)
        name = va["gnomad_sv_overlap_name"].split("&")[chosen_index]
        coords = va["gnomad_sv_overlap_coords"]
        if not coords:
            logging.warning("variant %s (%s) has no gnomad_sv_overlap_coords", va["variant_id"], name)
            results["no_coords"] += 1
            continue

        info = reader.get_info(name, coords.split("&")[chosen_index])
        if info is None:
            logging.warning("variant %s - %s not found in %s", va["variant_id"], name, reader.vcf_filename)
            results["record_not_found"] += 1
            continue

        # Proves the configured VCF is the one this version was annotated against - a different
        # gnomAD-SV release either names its records differently or reports a different frequency
        if float(info["AF"]) != lowest_af:
            logging.warning("variant %s - %s has AF %s in %s, annotation has %s",
                            va["variant_id"], name, info["AF"], reader.vcf_filename, lowest_af)
            results["af_mismatch"] += 1
            continue

        update = {}
        for source_field, columns in targets.items():
            value = info.get(source_field)
            for column in columns:
                update[column] = value

        if not dry_run:
            partition_qs.filter(variant_id=va["variant_id"]).update(**update)
        results["fixed"] += 1

        if results["checked"] % PROGRESS_INTERVAL == 0:
            logging.info("%s - checked %d, fixed %d", vav, results["checked"], results["fixed"])

    return results


class Command(BaseCommand):
    """ Re-pick the gnomAD-SV overlap record for SVs that overlapped more than one.

        The 'lowest_af' pick compared strings, so SVs whose overlapping records mixed scientific and
        decimal AF notation took their gnomAD columns from the most common record rather than the
        rarest. Reads the correct record back out of the configured gnomAD-SV VCF and rewrites them. """
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true",
                            help="Report what would change without writing")
        parser.add_argument("--variant-annotation-version", type=int,
                            help="Only process this VariantAnnotationVersion (default: all of them)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        qs = VariantAnnotationVersion.objects.all()
        if vav_id := options["variant_annotation_version"]:
            qs = qs.filter(pk=vav_id)

        totals = Counter()
        for vav in qs.order_by("pk"):
            logging.info("Processing %s (gnomAD-SV %s)", vav, vav.gnomad_sv or "version not recorded")
            with transaction.atomic():
                results = fix_variant_annotation_version(vav, dry_run=dry_run)
            logging.info("%s done: %s", vav, dict(results))
            totals.update(results)

        if dry_run:
            logging.info("Dry run - no changes written")
        logging.info("Totals: %s", dict(totals))
