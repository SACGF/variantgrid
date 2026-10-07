"""
The variant types page: how each kind of variant, by size, is stored, filtered in analyses, annotated,
lifted over and registered with ClinGen. Every cell is computed by the rule the code itself uses
(annotation.annotation_pipeline_routing, snpdb.bcftools_liftover, snpdb.clingen_allele), at each size
threshold in settings, so the page changes with them. Entry point: get_variant_type_rows.

Two of the thresholds are deliberately apart (#1358): a del/dup/inv is stored symbolic from
VARIANT_SYMBOLIC_ALT_SIZE, but annotated as a structural variant only from ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE.
"""
import itertools
from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings

from annotation.annotation_pipeline_routing import pipeline_type_for_alt
from annotation.models.models_enums import VariantAnnotationPipelineType
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.bcftools_liftover import (
    bcftools_liftover_skip_reason,
    bcftools_lifts_symbolic_as_explicit,
)
from snpdb.clingen_allele import ClinGenAlleleTooLargeException, clingen_check_variant_length
from snpdb.models import ClinGenAllele, VariantCoordinate
from snpdb.variant_filters import GENE_LEVEL_VARIANT_TYPES, VariantType, get_variant_type_label

_CHROM = "1"
_POSITION = 1_000_000


@dataclass(frozen=True)
class VariantKind:
    name: str
    # Coordinate as stored for a variant of this length - symbolic or explicit
    coordinate: Callable[[int], VariantCoordinate]
    sized: bool = True
    min_length: int = 1


def _symbolic(alt: str, svlen: int) -> VariantCoordinate:
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A", alt=alt, svlen=svlen)


def _stored_symbolic(min_length: int) -> Callable[[int], bool]:
    def _is_symbolic(length: int) -> bool:
        return settings.VARIANT_SYMBOLIC_ALT_ENABLED and length >= min_length
    return _is_symbolic


def _deletion(length: int) -> VariantCoordinate:
    # Babelfish makes a del symbolic when len(ref) (padding base included) > VARIANT_SYMBOLIC_ALT_SIZE
    if _stored_symbolic(settings.VARIANT_SYMBOLIC_ALT_SIZE)(length):
        return _symbolic(VCFSymbolicAllele.DEL, -length)
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A" * (length + 1), alt="A")


def _duplication(length: int) -> VariantCoordinate:
    if _stored_symbolic(settings.VARIANT_SYMBOLIC_ALT_SIZE)(length):
        return _symbolic(VCFSymbolicAllele.DUP, length)
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A", alt="A" * (length + 1))


def _inversion(length: int) -> VariantCoordinate:
    # An inv has no padding base, so is symbolic only when strictly longer (VariantCoordinate.from_vcf_coordinate)
    if _stored_symbolic(settings.VARIANT_SYMBOLIC_ALT_SIZE + 1)(length):
        return _symbolic(VCFSymbolicAllele.INV, length)
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A" * length, alt="T" * length)


def _insertion(length: int) -> VariantCoordinate:
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A", alt="A" + "C" * length)


def _complex_substitution(length: int) -> VariantCoordinate:
    return VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="AG" + "A" * length, alt="CT")


VARIANT_KINDS = [
    VariantKind("SNV", lambda _length: VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A", alt="T"),
                sized=False),
    VariantKind("Deletion", _deletion),
    VariantKind("Duplication", _duplication),
    VariantKind("Inversion", _inversion, min_length=2),  # 1 base is an SNV
    VariantKind("Insertion (not a duplication)", _insertion),
    VariantKind("Complex substitution", _complex_substitution),
    VariantKind("Copy number (<CNV>)", lambda length: _symbolic(VCFSymbolicAllele.CNV, length)),
    VariantKind("Insertion (<INS>)", lambda length: _symbolic(VCFSymbolicAllele.INS, length)),
]


def _size_thresholds() -> list[int]:
    """ Every length a rule below changes at. Each is also tried one longer, as some rules are > rather
        than >= - bands that turn out the same are merged """
    thresholds = {
        settings.VARIANT_SYMBOLIC_ALT_SIZE,
        settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE,
        settings.ANNOTATION_VEP_SV_MAX_SIZE,
        settings.LIFTOVER_BCFTOOLS_MAX_LENGTH,
        ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE,
        ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE // 2,  # dups count double
    }
    return sorted({length for t in thresholds if t for length in (t, t + 1) if length > 1})


def _length(vc: VariantCoordinate) -> int:
    if vc.is_symbolic:
        return abs(vc.svlen)
    return max(len(vc.ref), len(vc.alt)) - (len(vc.ref) != len(vc.alt))


def _stored_as(vc: VariantCoordinate) -> str:
    if vc.is_symbolic:
        return f"Symbolic {vc.alt} + SVLEN"
    return "Sequence (ref/alt)"


def _analysis_variant_type(vc: VariantCoordinate) -> str:
    """ The All Variants / analysis variant type (snpdb.variant_filters) """
    if vc.is_symbolic:
        return f"{get_variant_type_label(VariantType.SYMBOLIC)} ({get_variant_type_label(vc.alt)})"
    if len(vc.ref) == 1 and len(vc.alt) == 1:
        variant_type = VariantType.SNV
    elif len(vc.ref) == 1 or len(vc.alt) == 1:
        variant_type = VariantType.INDEL
    else:
        variant_type = VariantType.COMPLEX
    return get_variant_type_label(variant_type)


def _annotation(vc: VariantCoordinate) -> tuple[str, str]:
    pipeline_type = pipeline_type_for_alt(vc.alt, vc.svlen, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE)
    if pipeline_type == VariantAnnotationPipelineType.STANDARD:
        return pipeline_type.label, "gnomAD - exact match"
    if settings.ANNOTATION_VEP_SV_MAX_SIZE and _length(vc) > settings.ANNOTATION_VEP_SV_MAX_SIZE:
        return "Not annotated (over VEP's max SV size)", "-"
    same_type = ", same type" if settings.ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE else ""
    overlap = f"gnomAD-SV - overlapping SV (≥ {settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION:.0%}{same_type})"
    return pipeline_type.label, overlap


def _liftover(vc: VariantCoordinate) -> str:
    if not settings.LIFTOVER_BCFTOOLS_ENABLED:
        return "Disabled"
    if bcftools_lifts_symbolic_as_explicit(vc):
        return "Yes (as sequence)"
    if bcftools_liftover_skip_reason(vc):
        return "No"
    return "Yes"


def _clingen(vc: VariantCoordinate) -> str:
    """ As snpdb.clingen_allele._clingen_check_variant_coordinate_length """
    try:
        clingen_check_variant_length(str(vc), vc.max_sequence_length, is_dup=vc.alt == VCFSymbolicAllele.DUP)
    except ClinGenAlleleTooLargeException:
        return "No (too long)"
    return "Yes"


def _g_hgvs(vc: VariantCoordinate) -> str:
    if not vc.can_be_made_explicit:
        return "No"
    if vc.is_symbolic:
        return "Yes (from coordinates)"
    return "Yes"


def _cells(vc: VariantCoordinate) -> tuple:
    pipeline, population = _annotation(vc)
    return (_stored_as(vc), _analysis_variant_type(vc), pipeline, population,
            _liftover(vc), _clingen(vc), _g_hgvs(vc))


def _size_label(start: int, end: Optional[int]) -> str:
    if end is None:
        return f"≥ {start:,}bp"
    if end == start + 1:
        return f"{start:,}bp"
    return f"{start:,}-{end - 1:,}bp"


VARIANT_TYPE_COLUMNS = ["Stored as", "Analysis variant type", "Annotation pipeline", "Population frequency",
                        "Liftover (BCFtools)", "ClinGen Allele", "g.HGVS"]


def get_variant_type_rows() -> list[dict]:
    """ One row per kind and size band - adjacent bands that every column treats the same are merged """
    rows = []
    thresholds = _size_thresholds()
    for kind in VARIANT_KINDS:
        if not kind.sized:
            rows.append({"kind": kind.name, "size": "", "cells": _cells(kind.coordinate(1))})
            continue

        starts = [kind.min_length, *[t for t in thresholds if t > kind.min_length]]
        bands = [(start, _cells(kind.coordinate(start))) for start in starts]
        merged = [(next(group)[0], cells) for cells, group in itertools.groupby(bands, key=lambda band: band[1])]
        for (start, cells), next_band in itertools.zip_longest(merged, merged[1:]):
            end = next_band[0] if next_band else None
            rows.append({"kind": kind.name, "size": _size_label(start, end), "cells": cells})
    if settings.VARIANT_GENE_LEVEL_ENABLED:
        gene_level_types = ", ".join(get_variant_type_label(t) for t in GENE_LEVEL_VARIANT_TYPES)
        pipeline_type = pipeline_type_for_alt("", None, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE,
                                          is_gene_level=True)
        # A gene id where a coordinate goes, the same in every build (@see snpdb.gene_level_variants)
        rows.append({"kind": "Gene-level event", "size": "",
                     "cells": ("Gene id + event", gene_level_types, pipeline_type.label, "-",
                               "Not needed (same in every build)", "No", "No")})
    return rows
