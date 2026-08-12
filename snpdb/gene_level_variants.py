"""
Gene-level variants - why some Variants have no coordinate.

This is an unusual corner of the data model and not what anyone expects from `Variant`, so the
reasoning lives here rather than being spread across the places that have to know about it. If you
have arrived here from an `@see`, this file is the explanation.

Issue: https://github.com/SACGF/variantgrid/issues/1506
Design: claude/plans/fusion_variants_issue_1506_plan.md


## The problem

Some findings are real but have no position. A gene fusion is the case that forced this: a caller
reports "CD74-ROS1", and the biological identity is the *gene pair*, not a coordinate - the same
fusion turns up with different breakpoints in different reads of the same sample. (Coordinate-free
gene-level CNV, where a caller says only "JAK2 amplification", is the same shape and is designed for
but not yet built.)


## Why they are stored as Variants anyway

Deliberately a hack. The alternatives are worse:

  * A model of their own means writing every downstream feature a second time. Gene lists, filtering,
    grids, tags, classifications and analysis nodes are all built on `Variant`, and none of them
    would see a fusion.
  * A bare `Allele` (the original proposal on #1506) buys the classification chain and nothing else.
    Gene-list filtering and compound-het two-hit detection both join `VariantGeneOverlap`, which is a
    `Variant` FK, and sample-level presence *is* `CohortGenotype.variant`. With no `Variant` a fusion
    is not merely missing from a grid - it is unreachable from any analysis node.

So: accept one weird kind of `Variant`, and every existing feature works on fusions for free. That
held up in practice - gene lists, classifications and search all found fusions without learning what
a fusion is, and comp-het reaches them through the same `VariantGeneOverlap` rows a gene list uses.


## How it is arranged

  * One `Contig` (`SequenceRole.GENE_LEVEL`) shared by every genome build, using the constants below.
    Not the 5' gene's real chromosome: `AnalysisNode` narrows querysets by contig, an
    inter-chromosomal fusion names two chromosomes, and a gene list on the far partner would silently
    drop it. Sharing one contig across builds also means one `Variant` serves every build, so
    liftover has nothing to do.
  * `Locus.position` is a gene ID (a `genes.FusionGeneId` pk), NOT a coordinate. `Locus.ref` = 'N'.
  * The partner gene is encoded in the alt (`<FUSION:HGNC:nnn>`, @see `GeneLevelSymbolicAlt`), so identity
    hashes to its own `Sequence` row and the existing `(locus, alt, svlen)` unique constraint does the
    work with no new uniqueness machinery.
  * `svlen` is 0 rather than null, because Postgres treats nulls as distinct, so a null `svlen` would
    mean that constraint enforced nothing at all.


## What it costs, and where that is paid

These have no reference sequence, so anything that reads one has to be kept away from them.
`Variant.get_gene_level_q()` (and the `is_gene_level` instance twin) is the single predicate for
that - ClinGen registration, c.HGVS, and which annotation pipeline claims a variant all go through
it. Grep for it to find every such place.

VEP never sees a gene-level variant - it cannot parse the alt, and there is no coordinate to annotate
anyway. A third pipeline type (`VariantAnnotationPipelineType.GENE_LEVEL`) computes their annotation
from the gene identity instead, writing the `VariantGeneOverlap` rows for both partners that make
gene lists and comp-het work.

The VCF-writing paths need no check of their own: they build their contig list from
`GenomeBuild.standard_contigs`, which filters on `SequenceRole.ASSEMBLED_MOLECULE`, so the gene-level
contig is excluded before anything asks for a reference base. Worth keeping true - a writer that
switched to `contigs` rather than `standard_contigs` would start emitting them.
"""

# The contig every build shares. Deliberately non-genomic names so nothing mistakes it for a sequence,
# and a VG-namespaced accession so it cannot collide with a real one.
GENE_LEVEL_CONTIG_NAME = "GENE_LEVEL"
GENE_LEVEL_CONTIG_REFSEQ_ACCESSION = "VG_GENE_LEVEL"
# Well above any gene ID we will hand out - Variant.validate bounds position by the contig length
GENE_LEVEL_CONTIG_LENGTH = 1_000_000_000

# There is no reference base, so every gene-level Locus shares one
GENE_LEVEL_REF = "N"
# Not null - @see the module docstring
GENE_LEVEL_SVLEN = 0
