"""
Which genes of a GeneAnnotationRelease a genomic interval overlaps, resolved locally.

Two callers need this and neither can get it from VEP: long structural variants VEP skipped as
TOO_LONG (#1271, @see annotation.vcf_files.bulk_vep_vcf_annotation_inserter), and the breakpoints a
fusion caller writes, where position is the caller-spelling-independent evidence of which gene a
side is (@see genes.gene_fusions).

A TranscriptVersion keeps its extents inside the cdot JSON rather than in columns, so there is no
SQL predicate for "overlaps this range" - the transcripts of a contig go into an IntervalTree
instead. Trees are built per contig on first use: a release is ~180k transcripts and a fusion file
asks about a handful of positions, so building every contig up front would cost most of a gigabyte
to answer six lookups.
"""
import logging
import time
from dataclasses import dataclass
from typing import Optional

import intervaltree

from genes.models import GeneAnnotationRelease, TranscriptVersion
from snpdb.models import Contig, GenomeBuild, VariantCoordinate


@dataclass(frozen=True)
class GeneOverlap:
    """ One gene an interval falls in, as a GeneAnnotationRelease knows it. hgnc_id is what makes
        the same gene recognisable across releases and consortiums """
    gene_id: str
    symbol: Optional[str]
    hgnc_id: Optional[int]


class SVGeneOverlapResolver:
    """ The genes of one GeneAnnotationRelease overlapping a position or interval.

        Build one per release and keep it for the length of a run - the per-contig trees are the
        expensive part and are shared by every lookup that lands on that contig. """

    def __init__(self, gene_annotation_release: Optional[GeneAnnotationRelease]):
        self.gene_annotation_release = gene_annotation_release
        self._trees: dict[int, intervaltree.IntervalTree] = {}

    @classmethod
    def for_variant_annotation_version(cls, variant_annotation_version) -> 'SVGeneOverlapResolver':
        gene_annotation_release = variant_annotation_version.gene_annotation_release
        if gene_annotation_release is None:
            logging.warning("SVGeneOverlapResolver: no gene_annotation_release on %s",
                            variant_annotation_version)
        return cls(gene_annotation_release)

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        if self.gene_annotation_release:
            return self.gene_annotation_release.genome_build
        return None

    def _get_tree(self, contig: Contig) -> intervaltree.IntervalTree:
        if (tree := self._trees.get(contig.pk)) is not None:
            return tree

        start_time = time.monotonic()
        tree = intervaltree.IntervalTree()
        tv_qs = TranscriptVersion.objects.filter(
            releasetranscriptversion__release=self.gene_annotation_release,
            contig=contig,
        ).select_related("gene_version")

        count = 0
        for tv in tv_qs:
            try:
                start = tv.start
                end = tv.end
            except (KeyError, IndexError):
                continue
            if end <= start:
                # intervaltree treats zero-length intervals as empty
                end = start + 1
            gene_version = tv.gene_version
            tree.addi(start, end, GeneOverlap(gene_id=gene_version.gene_id,
                                              symbol=gene_version.gene_symbol_id,
                                              hgnc_id=gene_version.hgnc_id))
            count += 1

        self._trees[contig.pk] = tree
        logging.info("SVGeneOverlapResolver: %s %s - %d intervals in %.2fs",
                     self.gene_annotation_release, contig.name, count, time.monotonic() - start_time)
        return tree

    def get_gene_overlaps(self, chrom: str, start: int, end: Optional[int] = None) -> set[GeneOverlap]:
        """ chrom goes through GenomeBuild.chrom_contig_mappings, so 'chr3', '3' and NC_000003.11
            all reach the same contig - a fusion caller and a VCF spell a chromosome differently """

        if self.gene_annotation_release is None:
            return set()
        contig = self.genome_build.chrom_contig_mappings.get(chrom)
        if contig is None:
            return set()
        if end is None or end <= start:
            end = start + 1
        return {interval.data for interval in self._get_tree(contig).overlap(start, end)}

    def get_overlaps(self, variant_coordinate: VariantCoordinate) -> tuple[set[str], set[str]]:
        """ Returns (overlapping_symbols, overlapping_gene_ids) for the given variant. """
        symbols: set[str] = set()
        gene_ids: set[str] = set()
        for gene_overlap in self.get_gene_overlaps(variant_coordinate.chrom, variant_coordinate.position,
                                                   variant_coordinate.end):
            if gene_overlap.gene_id is not None:
                gene_ids.add(gene_overlap.gene_id)
            if gene_overlap.symbol is not None:
                symbols.add(gene_overlap.symbol)
        return symbols, gene_ids
