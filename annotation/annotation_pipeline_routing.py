"""
Which annotation pipeline a variant belongs to - the single source of truth, as a Q for querysets
(pipeline_type_variant_q) and per variant (pipeline_type_for_variant / pipeline_type_for_alt).

Storage and annotation have separate size cut-offs (#1358). A del/dup/inv is stored symbolic from
settings.VARIANT_SYMBOLIC_ALT_SIZE (50bp, the usual definition of a structural variant), but only goes to
the STRUCTURAL_VARIANT pipeline from settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE (1000bp). In between it
is annotated by the STANDARD pipeline, written to VEP as its explicit sequence
(annotation.annotation_run_files.write_qs_to_vcf). gnomAD's short-variant callset has exact allele
frequencies for indels into the hundreds of bases, where the SV pipeline only has gnomAD-SV overlap - and
the variants annotated before the storage cut-off came down keep their STANDARD annotation. <CNV>/<INS>
have no explicit form, so any size is an SV.
"""
from typing import Optional

from django.conf import settings
from django.db.models import Q

from annotation.models.models_enums import VariantAnnotationPipelineType
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models.models_variant import Sequence, Variant

# Symbolic alts that have an explicit sequence, so can be annotated as a small variant
EXPLICIT_SYMBOLIC_ALTS = (VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP, VCFSymbolicAllele.INV)


def symbolic_annotated_as_small(alt: str, svlen: Optional[int]) -> bool:
    """ A symbolic del/dup/inv short enough for the STANDARD pipeline. False for anything explicit """
    if alt not in EXPLICIT_SYMBOLIC_ALTS or svlen is None:
        return False
    return abs(svlen) < settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE


def pipeline_type_for_alt(alt: str, svlen: Optional[int], is_gene_level: bool = False) -> VariantAnnotationPipelineType:
    if is_gene_level:
        return VariantAnnotationPipelineType.GENE_LEVEL
    if Sequence.allele_is_symbolic(alt) and not symbolic_annotated_as_small(alt, svlen):
        return VariantAnnotationPipelineType.STRUCTURAL_VARIANT
    return VariantAnnotationPipelineType.STANDARD


def pipeline_type_for_variant(variant: Variant) -> VariantAnnotationPipelineType:
    return pipeline_type_for_alt(variant.alt.seq, variant.svlen, is_gene_level=variant.is_gene_level)


def _symbolic_annotated_as_small_q() -> Q:
    # svlen is negative for a VCF 4.3 <DEL>. alt_id IN (subquery) rather than a join to Sequence, which
    # costs whole-build scans their per-branch plans (@see snpdb.variant_filters._alt_in_q)
    min_size = settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE
    return Q(alt__in=Sequence.objects.filter(seq__in=EXPLICIT_SYMBOLIC_ALTS),
             svlen__gt=-min_size, svlen__lt=min_size)


def pipeline_type_variant_q(pipeline_type: VariantAnnotationPipelineType) -> Q:
    """ Which variants are subject to a given annotation pipeline type - the queryset twin of
        pipeline_type_for_alt. New pipeline types register their predicate here - a type's predicate may
        overlap another's. """

    q_sv = Variant.get_symbolic_q() & ~_symbolic_annotated_as_small_q()
    # Gene-level variants carry svlen=0 (so unique_together works - Postgres treats nulls as distinct),
    # which makes them look symbolic. VEP can't parse their alt and they have no coordinate anyway, so
    # subtract them from both VEP pipelines. @see snpdb.gene_level_variants
    q_gene_level = Variant.get_gene_level_q()
    if pipeline_type == VariantAnnotationPipelineType.STANDARD:
        return ~q_sv & ~q_gene_level
    elif pipeline_type in (VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                           VariantAnnotationPipelineType.ANNOTSV):
        # AnnotSV annotates the same variants VEP's SV pipeline does - a different tool, not a different
        # class of variant. The overlap the docstring above allows for.
        return q_sv & ~q_gene_level
    elif pipeline_type == VariantAnnotationPipelineType.GENE_LEVEL:
        return q_gene_level
    raise ValueError(f"Unrecognised {pipeline_type=}")
