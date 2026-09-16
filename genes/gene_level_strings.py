"""
Recognising a gene-level event written as one string, and resolving it.

A classification target arrives in the same field a c.HGVS does, so something has to decide which
path the value takes before either resolver runs: looks_gene_level answers that from the shape of
the string alone, so a value that names genes is never handed to the HGVS converter and never
reported as a bad transcript (@see classification.models.ImportedAlleleInfo).

resolve_gene_level_string then runs the three kinds' resolvers - fusion, whole-gene copy number,
splice - and hands back the identity, or the reason the kind that recognised the string refused it.
This is the only module that knows all three kinds; the resolvers themselves stay independent.
"""
from functools import partial
from typing import Optional

from genes.gene_copy_number import (
    COPY_NUMBER_STRING_PATTERN,
    resolve_gene_copy_number_string,
)
from genes.gene_fusions import FUSION_STRING_SEPARATOR, resolve_fusion_string
from genes.gene_level_resolver import GeneLevelResolution
from genes.gene_splice import SPLICE_STRING_PATTERN, resolve_splice_string
from snpdb.models import HGVS_UNCLEANED_PATTERN, GenomeBuild


def looks_gene_level(value: str) -> bool:
    """ Whether the value names genes ('BCR::ABL1', 'EGFR amplification', 'ARV7') rather than giving
        a coordinate. Shape only - it holds for a value nothing has resolved yet, and for one whose
        gene turned out to be a typo, which is what stops either being read as a broken HGVS """

    if not value or HGVS_UNCLEANED_PATTERN.search(value):
        return False
    return bool(FUSION_STRING_SEPARATOR.search(value)
                or COPY_NUMBER_STRING_PATTERN.match(value)
                or SPLICE_STRING_PATTERN.match(value))


def resolve_gene_level_string(value: str, genome_build: Optional[GenomeBuild] = None) -> GeneLevelResolution:
    """ The gene-level identity a written value names, or the reason it was refused. The build is
        the record's imported one - a junction named by its breakpoints means different junctions in
        different builds (@see genes.gene_splice). """

    if not looks_gene_level(value):
        return GeneLevelResolution.not_applicable()

    refusal: Optional[GeneLevelResolution] = None
    for resolve in (resolve_fusion_string, resolve_gene_copy_number_string,
                    partial(resolve_splice_string, genome_build=genome_build)):
        resolution = resolve(value)
        if resolution:
            return resolution
        if refusal is None and resolution.reason:
            refusal = resolution
    # a resolution is falsy until it resolved, so the first refusal is returned explicitly
    if refusal is not None:
        return refusal
    return GeneLevelResolution.refused(f"'{value}' names no gene-level event we can resolve")
