"""
The variant types page: how each kind of variant is stored and annotated, by size, and the size limits
of the features that have one. Every cell is computed by the rule the code itself uses
(annotation.annotation_pipeline_routing, snpdb.bcftools_liftover, snpdb.clingen_allele) at the thresholds
in settings, so the page changes with them. Entry points: get_variant_type_rows, get_size_limits,
get_unsupported.

Two of the thresholds are deliberately apart (#1358): a del/dup/inv is stored symbolic from
VARIANT_SYMBOLIC_ALT_SIZE, but annotated as a structural variant only from ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE.
Those two, and VEP's ceiling on a structural variant (ANNOTATION_VEP_SV_MAX_SIZE - above it a variant is
annotated with its overlapping genes only), are the sizes that split the table; a feature's own ceiling (liftover,
ClinGen) is a line under it, as it is the same for every kind.
"""
import itertools
from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings

from annotation.annotation_pipeline_routing import pipeline_type_for_alt
from annotation.models.models_enums import VariantAnnotationPipelineType
from library.genomics import format_bp
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import ClinGenAllele, VariantCoordinate

_CHROM = "1"
_POSITION = 1_000_000


@dataclass(frozen=True)
class VariantKind:
    name: str
    # Coordinate as stored for a variant of this length - symbolic or explicit
    coordinate: Callable[[int], VariantCoordinate]
    sized: bool = True
    min_length: int = 1
    symbolic_only: bool = False  # has no sequence form - stored only if its alt is in VARIANT_SYMBOLIC_ALT_VALID_TYPES


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


SNV = VariantKind("SNV", lambda _length: VariantCoordinate(chrom=_CHROM, position=_POSITION, ref="A", alt="T"),
                  sized=False)
DELETION = VariantKind("Deletion", _deletion)
DUPLICATION = VariantKind("Duplication", _duplication)
INVERSION = VariantKind("Inversion", _inversion, min_length=2)  # 1 base is an SNV
INSERTION = VariantKind("Insertion (not a duplication)", _insertion)
COMPLEX_SUBSTITUTION = VariantKind("Complex substitution", _complex_substitution)
CNV = VariantKind("Copy number (<CNV>)", lambda length: _symbolic(VCFSymbolicAllele.CNV, length), symbolic_only=True)
SYMBOLIC_INSERTION = VariantKind("Insertion (<INS>)", lambda length: _symbolic(VCFSymbolicAllele.INS, length),
                                 symbolic_only=True)

# Kinds that are usually treated alike share a row group - they get their own rows if settings make them differ
VARIANT_KIND_GROUPS = [
    ("SNV", [SNV]),
    ("Deletion, duplication, inversion", [DELETION, DUPLICATION, INVERSION]),
    ("Insertion, complex substitution", [INSERTION, COMPLEX_SUBSTITUTION]),
    ("Copy number <CNV>, insertion <INS>", [CNV, SYMBOLIC_INSERTION]),
]

# Symbolic alts left out of VARIANT_SYMBOLIC_ALT_VALID_TYPES for a reason worth giving
UNSUPPORTED_SYMBOLIC_REASONS = {
    VCFSymbolicAllele.INS: "Callers (eg Manta) often give no SVLEN and END=POS, so it would be stored with no "
                           "length and unrelated insertions at a position would become one variant.",
}


@dataclass(frozen=True)
class Outcome:
    """ What the table shows for a kind at a size - the key adjacent bands are merged on """
    symbolic: bool
    pipeline: VariantAnnotationPipelineType
    vep_too_long: bool = False

    @property
    def name(self) -> str:
        if not self.symbolic:
            return "Short"
        if self.pipeline == VariantAnnotationPipelineType.STANDARD:
            return "Small SV"
        if self.vep_too_long:
            return "Large SV"
        return "SV"


@dataclass(frozen=True)
class Band:
    start: int
    end: Optional[int]  # exclusive, None = open
    outcome: Outcome
    alts: tuple[str, ...]  # symbolic alts of the kinds in the band, for the Stored as cell
    end_inclusive: bool = False  # VEP's ceiling is the one limit a variant must be longer than, not reach


def _outcome(vc: VariantCoordinate) -> Outcome:
    pipeline_type = pipeline_type_for_alt(vc.alt, vc.svlen, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE)
    # As annotation_version_querysets.filter_vep_sv_max_size - VEP is never given these
    sv_max_size = settings.ANNOTATION_VEP_SV_MAX_SIZE
    vep_too_long = (pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT
                    and bool(sv_max_size) and vc.svlen is not None and abs(vc.svlen) > sv_max_size)
    return Outcome(symbolic=vc.is_symbolic, pipeline=pipeline_type, vep_too_long=vep_too_long)


def _size_thresholds() -> list[int]:
    """ The lengths where storage or annotation changes """
    thresholds = {settings.VARIANT_SYMBOLIC_ALT_SIZE, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE,
                  settings.ANNOTATION_VEP_SV_MAX_SIZE}
    return sorted(t for t in thresholds if t and t > 1)


def _kind_bands(kind: VariantKind) -> list[Band]:
    """ One band per threshold the kind's behaviour changes at. A band is probed one base above its start,
        as some rules are > rather than >= (an inversion has no padding base) - that off-by-one is a
        footnote, not a row """
    if not kind.sized:
        vc = kind.coordinate(1)
        return [Band(1, None, _outcome(vc), (vc.alt,) if vc.is_symbolic else ())]

    starts = [kind.min_length, *[t for t in _size_thresholds() if t > kind.min_length]]
    probes = [kind.min_length, *[start + 1 for start in starts[1:]]]
    coordinates = [(start, kind.coordinate(probe)) for start, probe in zip(starts, probes)]
    bands = []
    for outcome, group in itertools.groupby(coordinates, key=lambda sc: _outcome(sc[1])):
        start, vc = next(group)
        bands.append(Band(start, None, outcome, (vc.alt,) if vc.is_symbolic else ()))
    return [Band(band.start, next_band.start if next_band else None, band.outcome, band.alts,
                 end_inclusive=bool(next_band and next_band.outcome.vep_too_long))
            for band, next_band in itertools.zip_longest(bands, bands[1:])]


def _merge_bands(kinds_bands: list[list[Band]]) -> Optional[list[Band]]:
    """ The group's bands if every kind in it has the same ones (the first band's start may differ,
        as an inversion starts at 2), else None """
    first = kinds_bands[0]
    for bands in kinds_bands[1:]:
        if [(b.end, b.outcome) for b in bands] != [(b.end, b.outcome) for b in first]:
            return None
    merged = []
    for bands in zip(*kinds_bands):
        alts = tuple(dict.fromkeys(alt for band in bands for alt in band.alts))
        merged.append(Band(bands[0].start, bands[0].end, bands[0].outcome, alts, bands[0].end_inclusive))
    return merged


def _size_label(band: Band, sized: bool, first: bool, last: bool) -> tuple[str, str]:
    """ (name, range) eg ('Small SV', '50 bp to 1 kb') """
    if not sized:
        return "", ""
    if first and last:
        return "Any size", ""
    if first:
        return band.outcome.name, f"{'≤' if band.end_inclusive else '<'} {format_bp(band.end)}"
    if last:
        comparison = ">" if band.outcome.vep_too_long else "≥"
        return band.outcome.name, f"{comparison} {format_bp(band.start)}"
    return band.outcome.name, f"{format_bp(band.start)} to {format_bp(band.end)}"


def _annotation_cells(outcome: Outcome) -> tuple[str, str]:
    """ (pipeline, detail - population frequency, or what a variant VEP skips still gets) """
    if outcome.vep_too_long:
        detail = "Overlapping genes only - no consequence, transcripts, HGVS or gnomAD-SV"
        if settings.ANNOTATION_ANNOTSV_ENABLED:
            detail += " (AnnotSV still runs)"
        return f"{outcome.pipeline.label} - not annotated by VEP", detail
    if outcome.pipeline == VariantAnnotationPipelineType.STANDARD:
        pipeline = outcome.pipeline.label
        if outcome.symbolic:
            pipeline += " (given to VEP as its sequence)"
        return pipeline, "gnomAD - exact match"
    same_type = ", same type" if settings.ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE else ""
    overlap = f"gnomAD-SV - overlapping SV (≥ {settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION:.0%}{same_type})"
    return outcome.pipeline.label, overlap


def _band_row(band: Band, sized: bool, first: bool, last: bool) -> dict:
    name, size_range = _size_label(band, sized, first, last)
    stored_as = f"Symbolic {'/'.join(band.alts)} + SVLEN" if band.outcome.symbolic else "Sequence (ref/alt)"
    pipeline, detail = _annotation_cells(band.outcome)
    return {"name": name, "range": size_range, "stored_as": stored_as, "pipeline": pipeline, "detail": detail}


def _is_stored(kind: VariantKind) -> bool:
    """ As vcf_clean_alts - a symbolic-only kind is dropped at import unless its alt is a valid type """
    if not kind.symbolic_only:
        return True
    return settings.VARIANT_SYMBOLIC_ALT_ENABLED and kind.coordinate(1).alt in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES


def _kind_groups() -> list[tuple[str, list[VariantKind]]]:
    groups = []
    for name, kinds in VARIANT_KIND_GROUPS:
        stored = [k for k in kinds if _is_stored(k)]
        if stored:
            groups.append((name if stored == kinds else ", ".join(k.name for k in stored), stored))
    return groups


def get_variant_type_rows() -> list[dict]:
    """ One group per kind (or kinds that behave alike), with a row per size band - bands split only where
        storage or the annotation pipeline changes """
    groups = []
    for group_name, kinds in _kind_groups():
        kinds_bands = [_kind_bands(kind) for kind in kinds]
        merged = _merge_bands(kinds_bands)
        named = [(group_name, merged)] if merged else zip([k.name for k in kinds], kinds_bands)
        sized = kinds[0].sized
        for name, bands in named:
            rows = [_band_row(b, sized, i == 0, i == len(bands) - 1) for i, b in enumerate(bands)]
            groups.append({"kind": name, "bands": rows})
    if settings.VARIANT_GENE_LEVEL_ENABLED:
        pipeline_type = pipeline_type_for_alt("", None, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE,
                                              is_gene_level=True)
        # A gene id where a coordinate goes, the same in every build (@see snpdb.gene_level_variants)
        groups.append({"kind": "Gene-level event (fusion, copy number, splicing)",
                       "bands": [{"name": "", "range": "", "stored_as": "Gene id + event",
                                  "pipeline": pipeline_type.label, "detail": "-"}]})
    return groups


def _liftover_limit() -> str:
    if not settings.LIFTOVER_BCFTOOLS_ENABLED:
        return "Disabled on this server."
    max_length = settings.LIFTOVER_BCFTOOLS_MAX_LENGTH
    if not max_length:
        return "No size limit."
    text = (f"Up to {format_bp(max_length)}, by sequence - a symbolic variant within that is written out "
            f"in full for BCFtools.")
    if settings.LIFTOVER_BCFTOOLS_SYMBOLIC:
        return text + " Longer symbolic variants are lifted over symbolically; longer sequence variants are not."
    return text + " Longer variants are not lifted over."


def _clingen_limit() -> str:
    max_size = ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE
    # As Variant.clingen_allele_skip_reason: only kinds with a g.HGVS, and dups count double
    return (f"Up to {format_bp(max_size)}, a duplication up to {format_bp(max_size // 2)} as ClinGen counts "
            f"its sequence twice. {_no_g_hgvs_kinds()} are not registered.")


def _no_g_hgvs_kinds() -> str:
    alts = [kind.coordinate(1).alt for _, kinds in _kind_groups() for kind in kinds
            if not kind.coordinate(1).can_be_made_explicit]
    if settings.VARIANT_GENE_LEVEL_ENABLED:
        alts.append("gene-level events")
    if len(alts) > 1:
        return ", ".join(alts[:-1]) + " and " + alts[-1]
    return alts[0] if alts else ""


def _g_hgvs_limit() -> str:
    return (f"Any sequence variant; a symbolic {VCFSymbolicAllele.DEL}/{VCFSymbolicAllele.DUP}/"
            f"{VCFSymbolicAllele.INV} from its coordinates. {_no_g_hgvs_kinds()} have none.")


def get_size_limits() -> list[dict]:
    """ The ceilings that apply to every kind - what a feature can take, not what a variant is """
    limits = [
        ("Liftover (BCFtools)", _liftover_limit()),
        ("ClinGen Allele", _clingen_limit()),
        ("g.HGVS", _g_hgvs_limit()),
    ]
    return [{"name": name, "text": text} for name, text in limits if text]


def get_unsupported() -> list[dict]:
    """ What VCF import drops (vcf_clean_alts) - the table above is only what can be stored """
    if settings.VARIANT_SYMBOLIC_ALT_ENABLED:
        valid_types = settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES
        unsupported = [{"name": alt, "text": reason} for alt, reason in UNSUPPORTED_SYMBOLIC_REASONS.items()
                       if alt not in valid_types]
        unsupported.append({"name": "Other symbolic ALTs",
                            "text": "eg <BND>, <INS:ME> and <DEL:ME> - any symbolic ALT not in the table above. "
                                    "<DUP:TANDEM> is stored as <DUP>."})
    else:
        unsupported = [{"name": "Symbolic ALTs",
                        "text": "Symbolic variants are disabled on this server, so eg <DEL> and <CNV> are dropped."}]
    unsupported.append({"name": "Other bases",
                        "text": "eg N, IUPAC ambiguity codes or breakend notation - any ALT with a base other "
                                "than A, C, G or T."})
    return unsupported
