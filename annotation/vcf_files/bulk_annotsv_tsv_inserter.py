"""
AnnotSV TSV ingestion. Reads the AnnotSV-annotated TSV produced for a SV
AnnotationRun and updates the VariantAnnotation rows VEP already wrote.

Only "full" lines are persisted in this iteration. "split" lines (per-gene
overlaps) are counted for diagnostic logging only - persisting them needs a
gene-release mapping to our Gene table that AnnotSV's bundle-pinned snapshot
does not currently match. Tracked as a follow-up to issue #1533.
"""
import csv
import logging
import re
from typing import Any, Optional

from annotation.models.models import AnnotationRun, VariantAnnotation


# AnnotSV preserves the input VCF INFO column when run with -SVinputInfo 1.
# The dump VCF writes "variant_id=NNN" into INFO; we use that to join back.
_VARIANT_ID_RE = re.compile(r"(?:^|;)variant_id=(\d+)")


# AnnotSV TSV column name -> VariantAnnotation field name. Values verified
# against the bundle shipped with AnnotSV 3.5.8.
FULL_COLUMN_MAP: dict[str, str] = {
    "ACMG_class": "annotsv_acmg_class",
    "AnnotSV_ranking_score": "annotsv_acmg_score",
    "RE_gene": "annotsv_re_gene",
    "Repeat_type_left": "annotsv_repeat_type_left",
    "Repeat_type_right": "annotsv_repeat_type_right",
    "SegDup_left": "annotsv_segdup_left",
    "SegDup_right": "annotsv_segdup_right",
    "ENCODE_blacklist_left": "annotsv_encode_blacklist_left",
    "ENCODE_blacklist_right": "annotsv_encode_blacklist_right",
    "ENCODE_blacklist_characteristics_left": "annotsv_encode_blacklist_characteristics_left",
    "ENCODE_blacklist_characteristics_right": "annotsv_encode_blacklist_characteristics_right",
    "B_gain_AFmax": "annotsv_b_gain_af_max",
    "B_loss_AFmax": "annotsv_b_loss_af_max",
    "B_ins_AFmax": "annotsv_b_ins_af_max",
    "B_inv_AFmax": "annotsv_b_inv_af_max",
}

INT_FIELDS = {"annotsv_acmg_class"}
FLOAT_FIELDS = {
    "annotsv_acmg_score",
    "annotsv_b_gain_af_max",
    "annotsv_b_loss_af_max",
    "annotsv_b_ins_af_max",
    "annotsv_b_inv_af_max",
}

EMPTY_VALUES = {"", "NA", ".", "-1", "NaN", "nan"}


def _parse_value(field: str, raw: str) -> Optional[Any]:
    if raw is None:
        return None
    raw = raw.strip()
    if raw in EMPTY_VALUES:
        return None
    try:
        if field in INT_FIELDS:
            # AnnotSV ACMG_class is an int 1..5 but can occasionally be "NA".
            return int(raw)
        if field in FLOAT_FIELDS:
            return float(raw)
    except (TypeError, ValueError):
        return None
    return raw


def _extract_variant_id(row: dict[str, str]) -> Optional[int]:
    """ Read the variant PK from the preserved VCF INFO column. Falls back to
        the ID column if it happens to carry an integer (unit tests may write
        the variant id there directly). """
    info = row.get("INFO") or ""
    if m := _VARIANT_ID_RE.search(info):
        try:
            return int(m.group(1))
        except ValueError:
            return None
    raw = (row.get("ID") or "").strip()
    if raw and raw != ".":
        try:
            return int(raw)
        except ValueError:
            return None
    return None


def _row_to_update(row: dict[str, str]) -> dict[str, Any]:
    update: dict[str, Any] = {}
    for tsv_col, model_field in FULL_COLUMN_MAP.items():
        if tsv_col not in row:
            continue
        value = _parse_value(model_field, row[tsv_col])
        if value is not None:
            update[model_field] = value
    return update


def import_annotsv_tsv(annotation_run: AnnotationRun) -> int:
    """ Update existing VariantAnnotation rows for this run from AnnotSV's TSV.
        Returns count of rows updated. """
    if not annotation_run.annotsv_tsv_filename:
        return 0

    qs = VariantAnnotation.objects.filter(annotation_run=annotation_run)

    full_updates_by_variant: dict[int, dict[str, Any]] = {}
    split_count = 0
    full_count = 0
    skipped_no_id = 0

    with open(annotation_run.annotsv_tsv_filename, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mode = (row.get("Annotation_mode") or "").strip().lower()
            if mode == "split":
                split_count += 1
                continue
            if mode != "full":
                continue
            full_count += 1
            variant_id = _extract_variant_id(row)
            if variant_id is None:
                skipped_no_id += 1
                continue
            update = _row_to_update(row)
            if update:
                full_updates_by_variant[variant_id] = update

    updated = 0
    if full_updates_by_variant:
        existing = {va.variant_id: va for va in qs.filter(variant_id__in=full_updates_by_variant)}
        to_update: list[VariantAnnotation] = []
        fields_changed: set[str] = set()
        for variant_id, update in full_updates_by_variant.items():
            va = existing.get(variant_id)
            if va is None:
                continue
            for field, value in update.items():
                setattr(va, field, value)
                fields_changed.add(field)
            to_update.append(va)
        if to_update and fields_changed:
            VariantAnnotation.objects.bulk_update(to_update, sorted(fields_changed), batch_size=1000)
            updated = len(to_update)

    logging.info(
        "AnnotSV TSV import (run=%s): full=%d split=%d updated=%d no_id=%d",
        annotation_run.pk, full_count, split_count, updated, skipped_no_id,
    )

    annotation_run.annotsv_imported = True
    annotation_run.save()
    return updated
