"""
Which annotation pipeline a variant belongs to - the single source of truth, as a Q for querysets
(pipeline_type_variant_q) and per variant (pipeline_type_for_variant / pipeline_type_for_alt).

Kind first, then size. A gene-level variant is GENE_LEVEL before any size rule is read - it has no size, its
position is a gene id (@see snpdb.gene_level_variants). Everything below is about variants with a coordinate.

Storage and annotation have separate size cut-offs (#1358). A del/dup/inv is stored symbolic from
settings.VARIANT_SYMBOLIC_ALT_SIZE (50bp, the usual definition of a structural variant), but only goes to
the STRUCTURAL_VARIANT pipeline from VariantAnnotationVersion.structural_variant_min_size (1000bp, from
settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE). The cut-off is pinned on the version, so changing the setting makes
a new version rather than rerouting variants an existing one has already annotated. In between it is annotated
by the STANDARD pipeline, written to VEP as its explicit sequence
(annotation.annotation_run_files.write_qs_to_vcf). gnomAD's short-variant callset has exact allele
frequencies for indels into the hundreds of bases, where the SV pipeline only has gnomAD-SV overlap - and
the variants annotated before the storage cut-off came down keep their STANDARD annotation. <CNV>/<INS>
have no explicit form, so any size is an SV.

Symbolic storage depends on shape, so a long delins, an insertion that isn't a dup, or a substitution that
isn't an inversion is stored as its sequence however long it is (#1708). Those go to the STRUCTURAL_VARIANT
pipeline when their longer allele (ref or alt) reaches VariantAnnotationVersion.explicit_structural_variant_min_size,
written to VEP and AnnotSV as the sequence. The field is null on versions made before it existed, which route
explicit variants to STANDARD at any length - so nothing they annotated moves.
"""
from typing import Optional

from django.db.models import F, Q
from django.db.models.functions import Length

from annotation.enums import VariantAnnotationPipelineType
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models.models_variant import Sequence, Variant

# Symbolic alts that have an explicit sequence, so can be annotated as a small variant
EXPLICIT_SYMBOLIC_ALTS = (VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP, VCFSymbolicAllele.INV)


def symbolic_annotated_as_small(alt: str, svlen: Optional[int], structural_variant_min_size: int) -> bool:
    """ A symbolic del/dup/inv short enough for the STANDARD pipeline. False for anything explicit """
    if alt not in EXPLICIT_SYMBOLIC_ALTS or svlen is None:
        return False
    return abs(svlen) < structural_variant_min_size


def explicit_annotated_as_structural(ref: str, alt: str, explicit_structural_variant_min_size: Optional[int]) -> bool:
    """ An explicit (sequence) variant long enough for the STRUCTURAL_VARIANT pipeline. False for anything symbolic """
    if explicit_structural_variant_min_size is None or Sequence.allele_is_symbolic(alt):
        return False
    return max(len(ref), len(alt)) >= explicit_structural_variant_min_size


def pipeline_type_for_alt(ref: str, alt: str, svlen: Optional[int], structural_variant_min_size: int,
                          explicit_structural_variant_min_size: Optional[int]) -> VariantAnnotationPipelineType:
    """ STANDARD or STRUCTURAL_VARIANT for a coordinate variant. Gene-level variants are routed before this """
    if Sequence.allele_is_symbolic(alt):
        if not symbolic_annotated_as_small(alt, svlen, structural_variant_min_size):
            return VariantAnnotationPipelineType.STRUCTURAL_VARIANT
    elif explicit_annotated_as_structural(ref, alt, explicit_structural_variant_min_size):
        return VariantAnnotationPipelineType.STRUCTURAL_VARIANT
    return VariantAnnotationPipelineType.STANDARD


def pipeline_type_for_variant(variant: Variant, structural_variant_min_size: int,
                              explicit_structural_variant_min_size: Optional[int]) -> VariantAnnotationPipelineType:
    if variant.is_gene_level:
        return VariantAnnotationPipelineType.GENE_LEVEL
    return pipeline_type_for_alt(variant.locus.ref.seq, variant.alt.seq, variant.svlen, structural_variant_min_size,
                                 explicit_structural_variant_min_size)


def _symbolic_annotated_as_small_q(structural_variant_min_size: int) -> Q:
    # svlen is negative for a VCF 4.3 <DEL>. alt_id IN (subquery) rather than a join to Sequence, which
    # costs whole-build scans their per-branch plans (@see snpdb.variant_filters._alt_in_q)
    return Q(alt__in=Sequence.objects.filter(seq__in=EXPLICIT_SYMBOLIC_ALTS),
             svlen__gt=-structural_variant_min_size, svlen__lt=structural_variant_min_size)


def _explicit_annotated_as_structural_q(explicit_structural_variant_min_size: int) -> Q:
    # An explicit variant's end is position + len(ref) - 1, so the ref needs no join to Sequence
    q_long_ref = Q(end__gte=F("locus__position") + explicit_structural_variant_min_size - 1)
    long_alts = Sequence.objects.annotate(seq_length=Length("seq")).filter(
        seq_length__gte=explicit_structural_variant_min_size)
    return Q(svlen__isnull=True) & (q_long_ref | Q(alt__in=long_alts))


def pipeline_type_variant_q(pipeline_type: VariantAnnotationPipelineType, structural_variant_min_size: int,
                            explicit_structural_variant_min_size: Optional[int]) -> Q:
    """ Which variants are subject to a given annotation pipeline type - the queryset twin of
        pipeline_type_for_alt. New pipeline types register their predicate here - a type's predicate may
        overlap another's. """

    # A gene-level variant is GENE_LEVEL whatever its size, so it's routed before any size cut-off
    q_gene_level = Variant.get_gene_level_q()
    if pipeline_type == VariantAnnotationPipelineType.GENE_LEVEL:
        return q_gene_level

    q_sv = Variant.get_symbolic_q() & ~_symbolic_annotated_as_small_q(structural_variant_min_size)
    if explicit_structural_variant_min_size is not None:
        q_sv |= _explicit_annotated_as_structural_q(explicit_structural_variant_min_size)
    # Gene-level variants carry svlen=0 (so unique_together works - Postgres treats nulls as distinct),
    # which makes them look symbolic. VEP can't parse their alt and they have no coordinate anyway, so
    # subtract them from both VEP pipelines. @see snpdb.gene_level_variants
    if pipeline_type == VariantAnnotationPipelineType.STANDARD:
        return ~q_sv & ~q_gene_level
    elif pipeline_type in (VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                           VariantAnnotationPipelineType.ANNOTSV):
        # AnnotSV annotates the same variants VEP's SV pipeline does - a different tool, not a different
        # class of variant. The overlap the docstring above allows for.
        return q_sv & ~q_gene_level
    raise ValueError(f"Unrecognised {pipeline_type=}")
