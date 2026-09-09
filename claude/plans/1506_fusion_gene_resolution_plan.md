# Fusion gene identity: resolve by breakpoint, then HGNC previous/alias symbols

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: in progress

Follow-on to [#1506](https://github.com/SACGF/variantgrid/issues/1506) (gene fusions as Variants). The
storage design is `snpdb/gene_level_variants.py`; this plan changes only how each side of a fusion is
resolved to a gene and how the annotation finds release genes for it.

## The problem

Analysis 616 on vg-test2 holds ACPP::ETV1 (variant 16929494). A gene list node for ETV1 finds it, a gene
list node for ACP3 does not. ACPP is HGNC's previous symbol for ACP3 (HGNC:125).

The chain that breaks:

1. `genes/gene_matching.py:GeneSymbolMatcher.get_gene_symbol_id_and_alias_id` returns a direct `GeneSymbol`
   hit before consulting aliases. A `GeneSymbol` row `ACPP` exists (Ensembl still uses it), so the alias
   ACPP -> ACP3 is never applied.
2. `genes/gene_level_resolver.py:GeneLevelNameResolver.resolve_name` then finds no `HGNC` row with `gene_symbol=ACPP`,
   so `resolve_side` mints a custom `FusionGeneId` (pk 1000002, `gene_symbol=ACPP`, `hgnc=None`).
3. `annotation/gene_level_annotation.py:GeneLevelIdResolver` resolves the side through
   `GeneAnnotationRelease.genes_for_symbol("ACPP")`, which is empty in every RefSeq release (they know
   `ACP3` -> gene 55). Only the ETV1 side gets a `VariantGeneOverlap` row.
4. `analysis/models/nodes/filters/gene_list_node.py:GeneListNode` filters on `VariantGeneOverlap`, so ACP3
   (gene 55) matches nothing.

Of the 16 custom `FusionGeneId` rows on vg-test2, 11 are old symbols of real HGNC genes (NOTCH2NL, SEPT2,
ACPP, GAREM, FAM73A, ATP5G3, ZNF131, RPA3-AS1, RNVU1-11, C3orf83, MARCH9). Only the five clone-based
`RP11-…` names are genuinely without HGNC. The `FusionGeneId` docstring promises that a symbol which later
gains an HGNC entry has its `hgnc` filled in; nothing implements that.

The DRAGEN file carries breakpoints per call (ACPP::ETV1 is `chr3:132036420 -> chr7:13950932` and five other
pairs). They are kept only as text in `CohortGenotype.info` (`FUSION_OBS`). Position is the caller-spelling-
independent evidence of which gene a side is, so it becomes the primary resolution; the name stays as
tiebreaker and fallback.

## Data

```python
# genes/models/models_gene_fusion.py
class FusionGeneId(models.Model):
    symbol_str = TextField(unique=True, db_collation='case_insensitive')
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)
    # NEW: the release-independent genes this side is (Entrez id and ENSG for the same gene). Written by
    # breakpoint resolution at import. Annotation reads these first.
    genes = models.ManyToManyField(Gene, blank=True)
```

`GeneFusion`, `Variant`, `Locus` and the alt encoding are unchanged. Identity, once handed out, stays: a
custom pk that later gains an HGNC keeps its number and gets `hgnc` filled in, which is the policy the
docstring already states. New migration in `genes/`.

## Design

### 1. Name resolution learns HGNC previous and alias symbols

`genes/gene_level_resolver.py:GeneLevelNameResolver.resolve_name`, in order:

1. Current behaviour: matcher symbol -> `HGNC` row with that `gene_symbol` (approved preferred).
2. If no HGNC: the matcher's alias dict (`GeneSymbolAlias`), then a new upper-cased lookup built once per
   resolver from `HGNC.previous_symbols` and `HGNC.alias_symbols` (comma-separated text) -> HGNC, approved
   status preferred. A hit returns `(hgnc.gene_symbol_id, hgnc)`.
3. Otherwise as today: a known `GeneSymbol` with no HGNC, or nothing.

The lookup lives on the resolver as a `cached_property`, the same shape as `GeneSymbolMatcher._alias_dict`.

As implemented the order within step 2 is previous symbols, then `GeneSymbolAlias`, then alias symbols, and
the previous-symbol lookup also outranks a matcher hit that was itself an alias hop. Real data forces it:
`SEPT2` is `SEPTIN2`'s previous symbol *and* one of `SEPTIN6`'s aliases, `GeneSymbolAlias` holds a row for
each, and its upper-cased dict keeps only one - which on this box is `SEPTIN6`. A rename says two names are
one gene; an alias does not, so the rename decides. A name that is a current symbol in its own right is
still taken as written. All 11 old symbols in the problem statement now resolve, `SEPT2` to `SEPTIN2`.

### 2. Breakpoint resolution

`SVGeneOverlapResolver` is generalised to be built from a `GeneAnnotationRelease`; the existing
`VariantAnnotationVersion` entry point becomes a classmethod that reads the version's release, so
`fix_annotation_sv_overlaps` and the inserter keep their calls. Trees are keyed by `Contig` pk rather than
name and the lookup resolves a chromosome string through `GenomeBuild.chrom_contig_mappings`, so `chr3` and
`3` both hit.

As implemented it moved out of the annotation inserter to `genes/gene_overlaps.py:SVGeneOverlapResolver`:
`genes` cannot import `annotation.vcf_files.bulk_vep_vcf_annotation_inserter` at module level (that module
pulls in `upload`, and the genes/annotation cycle is real - `genes/CLAUDE.md`). The trees are also built per
contig on first use rather than per release up front: GRCh38 has six releases of ~177k transcripts, so
building whole releases to answer a dozen breakpoints would cost most of a gigabyte in a celery worker.

`GeneFusionResolver.resolve_side(cell, breakpoint=None, genome_build=None)`:

- `breakpoint` is the caller's `chr3:132036420` string. With a build, the resolver looks the position up in
  a tree per `GeneAnnotationRelease` of that build (built lazily, cached on the resolver, so one import pays
  once per release). Overlapping genes across all releases are grouped into candidate identities by HGNC
  (`GeneVersion.hgnc`), falling back to gene symbol for genes with none.
- The candidate group wins when the caller's name (after §1 resolution) names its HGNC or symbol. With no
  name match and exactly one group, that group wins. With several groups and no match, or no overlap at
  all, resolution falls back to the name alone (§1) and records no genes.
- The winning group gives the identity (its HGNC, or a custom id under its symbol) and every gene in the
  group is added to `FusionGeneId.genes`.

### 3. `FusionGeneId.get_or_create_for_symbol` is unchanged

Dropped along with the backfill. Nothing gives a custom row an `hgnc` after the fact any more, so a
"reuse the row that already has this HGNC" branch could never match - a re-load (§6) mints the
HGNC-numbered identity instead.

### 4. Annotation reads genes first, symbols second

`annotation/gene_level_annotation.py:GeneLevelIdResolver.get_release_gene_annotation`:

1. `fusion_gene_id.genes` restricted to genes in the release (`ReleaseGeneVersion`).
2. Else the symbol route as today, using `hgnc.gene_symbol_id` when `hgnc` is set and `gene_symbol_id`
   otherwise.

Transcript versions come from `GeneAnnotationRelease.transcript_versions_for_gene` over the resolved genes.
`overlapping_symbols` and `VariantAnnotation.symbol` use the approved symbol, so the grid shows ACP3::ETV1.

### 5. Import passes breakpoints and the build

`upload/tasks/import_dragen_tso500_all_fusions_task.py:_observations_by_variant_coordinate` runs before the
VCF exists, so the build is not yet resolved. Split the declared-or-source tail of
`upload/vcf/vcf_import.py:resolve_genome_build` into `resolve_genome_build_from_source(source, file_upload)`
and call it from the create-VCF step with the file's `# Source =` line. Pass `row.gene_a_breakpoint`,
`row.gene_b_breakpoint` and the build into `resolve_side`. With no build resolvable the import proceeds on
names, as it does today.

### 6. Re-load the caller's files

No backfill command. Identities minted before this change keep their custom pks, and what moves a
fusion onto its proper HGNC identity is re-loading the AllFusions.csv it came from: resolution now
mints the right one from the start.

Re-annotation needs no code either way. `VariantAnnotation`, `VariantTranscriptAnnotation` and
`VariantGeneOverlap` all hold `annotation_run` on `CASCADE`, so deleting the GENE_LEVEL
`AnnotationRun`s deletes their rows, and
`annotation/tasks/annotation_scheduler_task.py:_handle_variant_annotation_version` then recreates a
run for every range lock missing one - the same orphan scan that backfills a newly-enabled pipeline
(#720). A re-load doesn't even need that: its variants are new, so they get their own locks and runs.

What re-annotation alone cannot do is fix identity. A `FusionGeneId` with no `hgnc` and no `genes`
still resolves through a symbol no release carries, so re-running the pipeline over the existing rows
reproduces the same empty overlap. That is why the answer here is the re-load rather than a delete.

The old rows stay. Their variants keep whatever `CohortGenotype` they had, so a re-loaded sample's
fusions appear under the new identities and the superseded ones are reachable only from the old VCF.
Deleting that VCF leaves them without a sample, and so out of every analysis.

## Tests

- `genes/tests/test_gene_fusions.py`: ACPP resolves to ACP3 (HGNC previous symbol) when a `GeneSymbol` row
  `ACPP` exists and there is no `GeneSymbolAlias`; a breakpoint inside a transcript resolves to that gene
  regardless of the name written; a breakpoint in overlapping genes picks the one the name matches; a
  breakpoint in nothing falls back to the name; a rename outranks another gene's alias of the same name.
- `annotation/tests/test_gene_level_annotation.py`: a side with `genes` set writes `VariantGeneOverlap` for
  the release's member of that set and nothing for a symbol the release lacks; a gene list on the renamed
  side finds the fusion.
- `annotation/tests/test_sv_gene_overlaps.py`: keeps passing with the release-based constructor.
- `upload/tests/test_import_dragen_tso500_all_fusions.py`: the create-VCF step resolves the build from the
  source line and the written FUSION string carries the approved symbol.

## Out of scope

Merging two `FusionGeneId` rows (custom and HGNC) that describe one gene, and moving `Variant` rows between
identities. Clone-based names with no HGNC in any release stay custom ids with no genes.
