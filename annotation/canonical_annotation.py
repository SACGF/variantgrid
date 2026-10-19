"""
A variant's canonical annotation in the latest VariantAnnotationVersion of a build: get_canonical_transcript_annotation
(the canonical VariantTranscriptAnnotation) and get_canonical_c_hgvs (its c.HGVS with symbol, falling back to the
representative annotation for gene-level variants, which sit on no transcript).
"""
from typing import Optional

from annotation.models.models import VariantAnnotationVersion, VariantTranscriptAnnotation
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant


def get_canonical_transcript_annotation(variant: Variant, genome_build: GenomeBuild) -> Optional[VariantTranscriptAnnotation]:
    vav = VariantAnnotationVersion.latest(genome_build)
    return variant.varianttranscriptannotation_set.filter(version=vav, canonical=True).first()


def get_canonical_c_hgvs(variant: Variant, genome_build: GenomeBuild) -> Optional[str]:
    if cta := get_canonical_transcript_annotation(variant, genome_build):
        return cta.get_hgvs_c_with_symbol()
    if variant.is_gene_level:
        # The representative annotation carries the VICC gene-level nomenclature. @see snpdb.gene_level_variants
        vav = VariantAnnotationVersion.latest(genome_build)
        if va := variant.variantannotation_set.filter(version=vav).first():
            return va.get_hgvs_c_with_symbol()
    return None
