"""
Backfill individual VariantAnnotation columns from an annotated VCF (#1675).

Corrects or fills in named columns of an existing VariantAnnotationVersion without re-running VEP over
every variant. The driving case is #1673 - COSMIC keeps renaming the INFO field carrying its sample count
(CNT in v95/97, SAMPLE_COUNT in v99, GENOME_SCREEN_SAMPLE_COUNT in v101) so deployments hold cosmic_count
null, and re-annotating from scratch to recover one integer costs weeks.

Three steps, of which we own the first and third:

1. dump - write a VCF of the variants already annotated in a version, carrying variant_id in INFO
2. annotate - the operator annotates that VCF outside VariantGrid (`bcftools annotate` against a COSMIC
   VCF is the fast path; VEP works too)
3. import - read the annotated VCF back and UPDATE the named columns in batches

This is the mirror image of annotation.external_annotation at each step: it dumps variants that *are*
annotated, leaves AnnotationRun / range-lock state alone, UPDATEs named columns rather than INSERTing
whole rows, and its input deliberately comes from a different annotation source than the one that
produced the version.
"""
import csv
import gzip
import io
import logging
import os
import re
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Optional

import cyvcf2
from django.core.exceptions import FieldDoesNotExist
from django.db import connection, transaction
from django.db.models import BooleanField, Exists, OuterRef

from annotation import vep_columns as vep_columns_registry
from annotation.annotation_run_files import write_qs_to_vcf
from annotation.annotation_version_querysets import (
    get_queryset_for_annotation_version,
    get_variants_qs_for_annotation,
)
from annotation.enums import VariantAnnotationPipelineType
from annotation.models.models import VariantAnnotation, VariantTranscriptAnnotation
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_annotation import get_vep_skipped_variants_filename
from annotation.vep_field_formatters import EMPTY_VALUES, VEP_SEPARATOR
from annotation.vep_warnings import parse_skipped_variants_file
from genes.models.models_gene import TranscriptVersion
from library.django_utils import get_model_fields
from library.utils.file_utils import open_handle_gzip
from snpdb.models.models_variant import Variant
from snpdb.variants_to_vcf import VARIANT_ID
from upload.vcf.sql_copy_files import sql_copy_csv_file

# Records per COPY -> UPDATE cycle. A batch's rows are held in memory as CSV, so this trades memory for
# the number of times the partition is joined against.
DEFAULT_BATCH_SIZE = 1_000_000

# Temp table the batch is COPY'd into before being joined to the partitions. A row with a NULL
# transcript_version_id is the variant's own (INFO values and the PICK'd CSQ entry); one with a
# transcript_version_id is a CSQ entry's, for the transcript row the pipeline wrote for that Feature
TEMP_TABLE = "annotation_backfill_columns"
TEMP_TRANSCRIPT_COLUMN = "transcript_version_id"

# VEP writes its per-transcript annotations into this INFO field. A --csq target reads the PICK'd entry
# for the representative row and, when VariantTranscriptAnnotation holds the column, every entry for the
# transcript rows
CSQ_INFO_FIELD = "CSQ"
CSQ_FORMAT_SEPARATOR = "Format: "
CSQ_PICK_FIELD = "PICK"
CSQ_FEATURE_FIELD = "Feature"

IMPORT_PROGRESS_INTERVAL_SECONDS = 10

BACKFILL_FLAG_PREFIX = "backfilled_"


class BackfillColumnError(ValueError):
    """ A requested column can't be resolved against this version's VEPColumnDef registry entries. """


@dataclass(frozen=True)
class BackfillTarget:
    """ One VariantGrid column to backfill: where to read it from in the annotated VCF, how to clean the
        value, and which annotation tables hold it. """
    column: str  # VariantGridColumn id, which is also the model field / DB column name
    source_field: str  # INFO field (or CSQ field, when from_csq) to read it from
    from_csq: bool
    formatter: Optional[Callable]
    base_table_names: tuple[str, ...]

    @property
    def source_description(self) -> str:
        source = "CSQ" if self.from_csq else "INFO"
        return f"{source}/{self.source_field}"

    @property
    def in_transcript_table(self) -> bool:
        return VariantAnnotationVersion.TRANSCRIPT_ANNOTATION in self.base_table_names

    @property
    def per_transcript(self) -> bool:
        """ Read off every CSQ entry, each for the transcript row of its Feature """
        return self.from_csq and self.in_transcript_table


def _version_column_defs(variant_annotation_version: VariantAnnotationVersion) -> tuple:
    """ Registry entries that apply to this version. Pipeline type is left open - a variant-level column
        can be backfilled regardless of which pipeline originally wrote it. """
    return vep_columns_registry.filter_for(
        genome_build_name=variant_annotation_version.genome_build.name,
        columns_version=variant_annotation_version.columns_version,
        vep_version=variant_annotation_version.vep,
        cosmic_version=variant_annotation_version.cosmic,
        gnomad4_minor_version=variant_annotation_version.gnomad,
    )


def _base_table_names_for_column(column: str) -> tuple[str, ...]:
    """ Which annotation tables a column lives in: the representative transcript's, plus the per-transcript
        table when VariantTranscriptAnnotation holds the column. An INFO field is variant-level - the same
        value the pipeline copies into every transcript row; a CSQ field is read per entry, each written to
        the transcript row of its Feature. """
    if column not in set(get_model_fields(VariantAnnotation)):
        raise BackfillColumnError(f"Column '{column}' is not a field of VariantAnnotation")

    base_table_names = [VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION]
    if column in set(get_model_fields(VariantTranscriptAnnotation)):
        base_table_names.append(VariantAnnotationVersion.TRANSCRIPT_ANNOTATION)
    return tuple(base_table_names)


def resolve_backfill_columns(variant_annotation_version: VariantAnnotationVersion,
                             column_ids: Iterable[str],
                             info_overrides: Optional[dict[str, str]] = None,
                             csq_overrides: Optional[dict[str, str]] = None) -> list[BackfillTarget]:
    """ Turn requested VariantGrid column ids into BackfillTargets, using the VEPColumnDef registry for
        the formatter and the default source field. Anything that can't be resolved raises, so a typo or
        a column that doesn't exist at this version's columns_version fails before anything is dumped. """
    info_overrides = dict(info_overrides or {})
    csq_overrides = dict(csq_overrides or {})
    column_ids = list(column_ids)

    if both := set(info_overrides) & set(csq_overrides):
        raise BackfillColumnError(f"Columns given both --info and --csq sources: {', '.join(sorted(both))}")
    if unknown := (set(info_overrides) | set(csq_overrides)) - set(column_ids):
        raise BackfillColumnError(f"--info/--csq given for columns not in --columns: {', '.join(sorted(unknown))}")

    version_defs = _version_column_defs(variant_annotation_version)
    targets = []
    for column in column_ids:
        column_defs = [c for c in version_defs if column in c.variant_grid_columns]
        if not column_defs:
            known = [c for c in vep_columns_registry.VEP_COLUMNS if column in c.variant_grid_columns]
            if known:
                reason = (f"it is not populated by {variant_annotation_version} "
                          f"(columns_version={variant_annotation_version.columns_version}, "
                          f"vep={variant_annotation_version.vep})")
            else:
                reason = "no VEPColumnDef writes it"
            raise BackfillColumnError(f"Can't backfill column '{column}' - {reason}")

        formatters = {c.formatter for c in column_defs}
        if len(formatters) > 1:
            raise BackfillColumnError(f"Column '{column}' has {len(formatters)} different registry "
                                      f"formatters for {variant_annotation_version} - can't tell which "
                                      f"one wrote the existing values")
        formatter = formatters.pop()

        if column in csq_overrides:
            source_field, from_csq = csq_overrides[column], True
        elif column in info_overrides:
            source_field, from_csq = info_overrides[column], False
        else:
            source_fields = {c.vep_info_field for c in column_defs if c.vep_info_field}
            if not source_fields:
                raise BackfillColumnError(f"Column '{column}' is calculated at insert time, not read from "
                                          f"a VEP field - name its source with --info or --csq")
            if len(source_fields) > 1:
                fields = ", ".join(sorted(source_fields))
                raise BackfillColumnError(f"Column '{column}' is written from several fields ({fields}) - "
                                          f"name the one to read with --info or --csq")
            source_field, from_csq = source_fields.pop(), False

        targets.append(BackfillTarget(column=column, source_field=source_field, from_csq=from_csq,
                                      formatter=formatter,
                                      base_table_names=_base_table_names_for_column(column)))
    return targets


def _validate_annotation_columns(columns: Iterable[str]) -> list[str]:
    columns = list(columns or [])
    if unknown := set(columns) - set(get_model_fields(VariantAnnotation)):
        raise BackfillColumnError(f"Not fields of VariantAnnotation: {', '.join(sorted(unknown))}")
    return columns


def dump_annotated_variants(variant_annotation_version: VariantAnnotationVersion, output_filename: str,
                            targets: Optional[list[BackfillTarget]] = None, only_missing: bool = False,
                            not_null_columns: Optional[Iterable[str]] = None, consequence: Optional[str] = None,
                            min_variant_id: Optional[int] = None, max_variant_id: Optional[int] = None,
                            pipeline_type: Optional[VariantAnnotationPipelineType] = None) -> int:
    """ Write the variants annotated in this version to a VCF carrying variant_id in INFO, ready to be
        annotated externally. Returns the number of records written.

        only_missing restricts to rows where every target column is null - the common backfill shape, and
        much smaller than a full dump.

        not_null_columns restricts to rows that already have a value for those columns. A source writes
        all of its columns together, so a sibling column being null rules the variant out - #1822's Open
        Targets score only needed the 465k variants the plugin matched, not the version's 8.4M.

        consequence restricts to variants with that SO term on any transcript row - for a VEP plugin that
        only annotates some consequences, the variants whose transcript rows (and representative row) a
        --csq import can fill (#2426).

        pipeline_type restricts to the variants that pipeline annotates - STANDARD for a VEP plugin column,
        so the structural variants the STANDARD VEP command never sees stay out. Gene-level variants have
        annotation but no coordinate, so are never dumped. """
    annotation_version = variant_annotation_version.get_any_annotation_version()
    qs = get_variants_qs_for_annotation(annotation_version, pipeline_type=pipeline_type, annotated=True,
                                        min_variant_id=min_variant_id, max_variant_id=max_variant_id)
    qs = qs.exclude(Variant.get_gene_level_q())
    # get_variants_qs_for_annotation(annotated=True) only stops it excluding annotated variants - to dump
    # the ones that *have* annotation for this version we have to join to the partition ourselves.
    # One filter() call, so the multi-valued relation is joined once (a column's __isnull=True on its own
    # would also match variants with no annotation row at all)
    annotation_filter = {"variantannotation__isnull": False}
    if only_missing:
        if not targets:
            raise ValueError("only_missing requires the resolved targets")
        for target in targets:
            annotation_filter[f"variantannotation__{target.column}__isnull"] = True
    for column in _validate_annotation_columns(not_null_columns):
        annotation_filter[f"variantannotation__{column}__isnull"] = False
    qs = qs.filter(**annotation_filter)
    if consequence:
        # consequence is VEP's '&'-joined SO terms. An EXISTS rather than a join, so a variant with the
        # term on several transcripts is dumped once
        consequence_regex = rf"(^|{VEP_SEPARATOR}){re.escape(consequence)}({VEP_SEPARATOR}|$)"
        transcript_rows = get_queryset_for_annotation_version(VariantTranscriptAnnotation, annotation_version) \
            .filter(variant_id=OuterRef("pk"), consequence__regex=consequence_regex)
        qs = qs.filter(Exists(transcript_rows))

    # Matched by sequence against the source's file, like the STANDARD dump - an SV stays symbolic
    return write_qs_to_vcf(output_filename, variant_annotation_version.genome_build, qs,
                           structural_variant_min_size=variant_annotation_version.structural_variant_min_size)


def _count_vcf_records(vcf_filename: str) -> int:
    try:
        with open_handle_gzip(vcf_filename, "rt") as f:
            return sum(1 for line in f if not line.startswith("#"))
    except (EOFError, gzip.BadGzipFile) as e:
        raise BackfillColumnError(f"'{vcf_filename}' is truncated or corrupt: {e}") from e


def check_annotated_vcf_complete(dump_vcf_filename: str, annotated_vcf_filename: str) -> int:
    """ Every dumped record must come back annotated or listed in VEP's --skipped_variants_file: VEP can exit 0
        having written a short VCF (#1701), and a backfill flagged done over one leaves the rest null for good """
    dump_count = _count_vcf_records(dump_vcf_filename)
    annotated_count = _count_vcf_records(annotated_vcf_filename)
    skipped_count = 0
    skipped_filename = get_vep_skipped_variants_filename(annotated_vcf_filename)
    if os.path.exists(skipped_filename):
        skipped_count = len(parse_skipped_variants_file(skipped_filename))
    if annotated_count + skipped_count != dump_count:
        raise BackfillColumnError(f"'{annotated_vcf_filename}' has {annotated_count} records and VEP skipped "
                                  f"{skipped_count}, but '{dump_vcf_filename}' has {dump_count} - the annotation "
                                  f"stopped early. Re-run it")
    return dump_count


def validate_backfill_flag(flag: str):
    """ A VariantAnnotationVersion backfilled_* flag - what readers check to tell an un-backfilled null
        from a real one """
    if not flag.startswith(BACKFILL_FLAG_PREFIX):
        raise BackfillColumnError(f"'{flag}' is not a {BACKFILL_FLAG_PREFIX}* flag")
    try:
        field = VariantAnnotationVersion._meta.get_field(flag)
    except FieldDoesNotExist:
        raise BackfillColumnError(f"VariantAnnotationVersion has no field '{flag}'") from None
    if not isinstance(field, BooleanField):
        raise BackfillColumnError(f"VariantAnnotationVersion.{flag} is not a flag")


def mark_backfilled(variant_annotation_version: VariantAnnotationVersion, flag: str):
    validate_backfill_flag(flag)
    VariantAnnotationVersion.objects.filter(pk=variant_annotation_version.pk).update(**{flag: True})


class _CSQReader:
    """ Reads the CSQ entries of a VEP-annotated record: the PICK'd one for the representative row, and
        each entry's Feature for the transcript rows. """

    def __init__(self, reader: cyvcf2.VCF, per_transcript: bool):
        try:
            description = reader.get_header_type(CSQ_INFO_FIELD)["Description"]
        except KeyError:
            raise BackfillColumnError(f"VCF has no {CSQ_INFO_FIELD} INFO field - it wasn't annotated by "
                                      f"VEP, so --csq has nothing to read") from None
        _, _, csq_format = description.partition(CSQ_FORMAT_SEPARATOR)
        if not csq_format:
            raise BackfillColumnError(f"VCF {CSQ_INFO_FIELD} header has no '{CSQ_FORMAT_SEPARATOR}' - "
                                      f"can't tell which field is which")
        self.fields = csq_format.strip().strip('"').split("|")
        self.pick_index = self.fields.index(CSQ_PICK_FIELD) if CSQ_PICK_FIELD in self.fields else None
        self.feature_index = None
        if per_transcript:
            if CSQ_FEATURE_FIELD not in self.fields:
                raise BackfillColumnError(f"VCF {CSQ_INFO_FIELD} has no {CSQ_FEATURE_FIELD} field, so its entries "
                                          f"can't be matched to transcript rows - run VEP with "
                                          f"'--fields {CSQ_FEATURE_FIELD},...'")
            self.feature_index = self.fields.index(CSQ_FEATURE_FIELD)

    def index_of(self, source_field: str) -> int:
        try:
            return self.fields.index(source_field)
        except ValueError:
            raise BackfillColumnError(f"'{source_field}' is not a {CSQ_INFO_FIELD} field of this VCF "
                                      f"({', '.join(self.fields)})") from None

    def entries(self, variant) -> list[list[str]]:
        csq = variant.INFO.get(CSQ_INFO_FIELD)
        if not csq:
            return []
        return [entry.split("|") for entry in csq.split(",")]

    def picked_entry(self, entries: list[list[str]]) -> Optional[list[str]]:
        if not entries:
            return None
        if self.pick_index is not None:
            for entry in entries:
                if entry[self.pick_index] == "1":
                    return entry
        return entries[0]

    def feature(self, entry: list[str]) -> Optional[str]:
        if self.feature_index is None or self.feature_index >= len(entry):
            return None
        return entry[self.feature_index]


class _TranscriptVersionResolver:
    """ VEP's Feature -> TranscriptVersion pk, as the pipeline inserter links a transcript row: a versioned
        accession (the pipeline passes --transcript_version, so Ensembl ids carry one like RefSeq's) looked
        up among the version's build and consortium. A Feature the database has no version for matches
        no row, as the inserter left that row unlinked. """

    def __init__(self, variant_annotation_version: VariantAnnotationVersion):
        self.transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id(
            variant_annotation_version.genome_build, variant_annotation_version.annotation_consortium)

    def resolve(self, feature: Optional[str]) -> Optional[int]:
        if not feature:
            return None
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(feature)
        return self.transcript_versions_by_id.get(transcript_id, {}).get(version)


def normalise_source_value(value) -> Optional[str]:
    """ cyvcf2 hands back an int/float/str for a single INFO value and a tuple for a multi-valued one,
        while the registry formatters take the '&'-joined string VEP produces. Normalising to that string
        is what lets one code path serve both bcftools and VEP input. """
    if value is None or value is False:
        return None
    if value is True:  # INFO flag
        return "1"
    if isinstance(value, (tuple, list)):
        values = [str(v) for v in value if v is not None and str(v) not in EMPTY_VALUES]
        value = VEP_SEPARATOR.join(values)
    else:
        value = str(value)
    value = value.strip()
    if value in EMPTY_VALUES:
        return None
    return value


def _format_value(target: BackfillTarget, raw_value) -> Optional[str]:
    value = normalise_source_value(raw_value)
    if value is None:
        return None
    if target.formatter:
        value = target.formatter(value)
    return value


def _csv_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _create_temp_table(variant_annotation_version: VariantAnnotationVersion, targets: list[BackfillTarget]):
    """ Temp table shaped like the columns we're about to write, taking its types from the version's own
        partition so a COPY of formatted values lands in exactly the type the UPDATE compares against. """
    partition = variant_annotation_version.get_partition_table(
        base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
    columns = ", ".join(f'"{t.column}"' for t in targets)
    with connection.cursor() as cursor:
        cursor.execute(f"DROP TABLE IF EXISTS {TEMP_TABLE}")
        cursor.execute(f"CREATE TEMP TABLE {TEMP_TABLE} AS "
                       f"SELECT variant_id, {TEMP_TRANSCRIPT_COLUMN}, {columns} FROM {partition} WITH NO DATA")


def _drop_temp_table():
    with connection.cursor() as cursor:
        cursor.execute(f"DROP TABLE IF EXISTS {TEMP_TABLE}")


def _new_value_sql(column: str, clear_missing: bool) -> str:
    """ What the column becomes. By default a source with no value leaves the existing value alone, so a
        partial annotation file fills in what it knows; clear_missing makes the file authoritative for
        the dumped variants, writing NULL where the source is silent. """
    if clear_missing:
        return f't."{column}"'
    return f'COALESCE(t."{column}", a."{column}")'


def _batch_updates(targets: list[BackfillTarget]) -> list[tuple[str, str, list[str]]]:
    """ (base table, join from the partition `a` to the temp table `t`, columns) for each UPDATE a batch
        runs: the representative row and a variant-level column's transcript rows take the variant's own
        temp row, a per-transcript column's transcript rows each take their Feature's """
    variant_row = f"a.variant_id = t.variant_id AND t.{TEMP_TRANSCRIPT_COLUMN} IS NULL"
    transcript_row = f"a.variant_id = t.variant_id AND a.{TEMP_TRANSCRIPT_COLUMN} = t.{TEMP_TRANSCRIPT_COLUMN}"
    transcript_table = VariantAnnotationVersion.TRANSCRIPT_ANNOTATION
    updates = [
        (VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION, variant_row, [t.column for t in targets]),
        (transcript_table, variant_row, [t.column for t in targets if t.in_transcript_table and not t.per_transcript]),
        (transcript_table, transcript_row, [t.column for t in targets if t.per_transcript]),
    ]
    return [update for update in updates if update[2]]


def _update_batch(variant_annotation_version: VariantAnnotationVersion, targets: list[BackfillTarget],
                  clear_missing: bool, dry_run: bool) -> int:
    """ Join the temp table to each partition holding the targets. Returns rows changed (or, for dry_run,
        rows that would change). """
    rows_changed = 0
    with connection.cursor() as cursor:
        cursor.execute(f"ANALYZE {TEMP_TABLE}")
        for base_table_name, join, columns in _batch_updates(targets):
            partition = variant_annotation_version.get_partition_table(base_table_name=base_table_name)
            # IS DISTINCT FROM lets postgres skip rewriting rows whose value is unchanged - for COSMIC
            # only a small fraction of variants carry a count, which is most of why this is quick
            changed = " OR ".join(f'a."{c}" IS DISTINCT FROM {_new_value_sql(c, clear_missing)}'
                                  for c in columns)
            if dry_run:
                cursor.execute(f"SELECT count(*) FROM {partition} a "
                               f"INNER JOIN {TEMP_TABLE} t ON {join} "
                               f"WHERE {changed}")
                rows_changed += cursor.fetchone()[0]
            else:
                set_clause = ", ".join(f'"{c}" = {_new_value_sql(c, clear_missing)}' for c in columns)
                cursor.execute(f"UPDATE {partition} a SET {set_clause} FROM {TEMP_TABLE} t "
                               f"WHERE {join} AND ({changed})")
                rows_changed += cursor.rowcount
        cursor.execute(f"TRUNCATE {TEMP_TABLE}")
    return rows_changed


def _copy_and_update(variant_annotation_version: VariantAnnotationVersion, targets: list[BackfillTarget],
                     buffer: io.StringIO, clear_missing: bool, dry_run: bool) -> int:
    columns = ["variant_id", TEMP_TRANSCRIPT_COLUMN] + [t.column for t in targets]
    data = io.BytesIO(buffer.getvalue().encode())
    with transaction.atomic():
        sql_copy_csv_file(data, TEMP_TABLE, columns, quote='"')
        rows_changed = _update_batch(variant_annotation_version, targets, clear_missing, dry_run)
        if dry_run:
            transaction.set_rollback(True)
    return rows_changed


def import_backfill_vcf(variant_annotation_version: VariantAnnotationVersion, vcf_filename: str,
                        targets: list[BackfillTarget], batch_size: int = DEFAULT_BATCH_SIZE,
                        clear_missing: bool = False, dry_run: bool = False,
                        emit: Optional[Callable] = None) -> dict:
    """ Read an annotated VCF and update the target columns of the version's partitions in batches.
        Returns counts of records read / rows COPY'd / rows changed. """
    def _emit(message):
        if emit:
            emit(message)
        else:
            logging.info(message)

    reader = cyvcf2.VCF(vcf_filename)
    csq_reader = None
    csq_indexes = {}
    transcript_resolver = None
    per_transcript_targets = [t for t in targets if t.per_transcript]
    if any(t.from_csq for t in targets):
        csq_reader = _CSQReader(reader, per_transcript=bool(per_transcript_targets))
        csq_indexes = {t.column: csq_reader.index_of(t.source_field) for t in targets if t.from_csq}
    if per_transcript_targets:
        transcript_resolver = _TranscriptVersionResolver(variant_annotation_version)

    def _entry_value(entry: Optional[list[str]], target: BackfillTarget) -> Optional[str]:
        index = csq_indexes[target.column]
        return _format_value(target, entry[index] if entry and index < len(entry) else None)

    records_read = 0
    rows_copied = 0
    transcript_rows_copied = 0
    rows_changed = 0
    batch_rows = 0
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    last_progress = time.time()

    _create_temp_table(variant_annotation_version, targets)
    try:
        for variant in reader:
            records_read += 1
            variant_id = variant.INFO.get(VARIANT_ID)
            if variant_id is None:
                raise BackfillColumnError(f"{vcf_filename} record {variant.CHROM}:{variant.POS} has no "
                                          f"'{VARIANT_ID}' INFO field - it wasn't dumped by "
                                          f"annotation_backfill_columns --dump")

            csq_entries = csq_reader.entries(variant) if csq_reader else []
            picked_entry = csq_reader.picked_entry(csq_entries) if csq_reader else None
            values = []
            for target in targets:
                if target.from_csq:
                    values.append(_entry_value(picked_entry, target))
                else:
                    values.append(_format_value(target, variant.INFO.get(target.source_field)))

            # Without clear_missing a record the source says nothing about can't change anything, so it's
            # not worth COPYing - which for COSMIC is the great majority of the dump
            if clear_missing or any(v is not None for v in values):
                writer.writerow([variant_id, ""] + [_csv_value(v) for v in values])
                batch_rows += 1

            # One row per CSQ entry whose Feature is a transcript row, carrying only the per-transcript
            # columns. A plugin that annotates few transcripts (UTRAnnotator's 5' UTRs) leaves most empty
            for entry in csq_entries:
                transcript_version_id = transcript_resolver.resolve(csq_reader.feature(entry)) \
                    if transcript_resolver else None
                if transcript_version_id is None:
                    continue
                entry_values = {t.column: _entry_value(entry, t) for t in per_transcript_targets}
                if clear_missing or any(v is not None for v in entry_values.values()):
                    writer.writerow([variant_id, transcript_version_id]
                                    + [_csv_value(entry_values.get(t.column)) for t in targets])
                    batch_rows += 1
                    transcript_rows_copied += 1

            if batch_rows >= batch_size:
                rows_changed += _copy_and_update(variant_annotation_version, targets, buffer,
                                                 clear_missing, dry_run)
                rows_copied += batch_rows
                _emit(f"Read {records_read} records, {rows_copied} written, {rows_changed} rows "
                      f"{'would change' if dry_run else 'changed'}")
                buffer = io.StringIO()
                writer = csv.writer(buffer)
                batch_rows = 0
                last_progress = time.time()
            elif time.time() - last_progress > IMPORT_PROGRESS_INTERVAL_SECONDS:
                _emit(f"Read {records_read} records ({batch_rows} pending in batch)")
                last_progress = time.time()

        if batch_rows:
            rows_changed += _copy_and_update(variant_annotation_version, targets, buffer,
                                             clear_missing, dry_run)
            rows_copied += batch_rows
    finally:
        reader.close()
        _drop_temp_table()

    columns = ", ".join(t.column for t in targets)
    transcript_rows = f" ({transcript_rows_copied} for transcript rows)" if per_transcript_targets else ""
    _emit(f"Backfilled {columns}: read {records_read} records, wrote {rows_copied}{transcript_rows}, "
          f"{rows_changed} rows {'would change' if dry_run else 'changed'}")
    return {"records_read": records_read, "rows_copied": rows_copied,
            "transcript_rows_copied": transcript_rows_copied, "rows_changed": rows_changed}
