# Gene fusions as variants — #1506 phase 1

Design for [#1506](https://github.com/SACGF/variantgrid/issues/1506) phase 1: storing gene fusions (and,
later, coordinate-free gene-level CNV) as `Variant` rows, and the `AllFusions.csv` parser that produces
them. This is Phase 5 of [`tso500_overall_plan.md`](tso500_overall_plan.md); test data and its gotchas
are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).

Fusion equivalence and discordance — the part the Shariant user group called research-level — stays out
of scope. Somatic classifications already land on `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED`, so
nothing regresses by leaving it.

---

## The premise corrections this plan is built on

Three things changed between the design in #1506's comments and this plan.

**Fusions get a `Variant`, not a bare `Allele`.** The earliest design on the issue (April) gave a
`GeneFusion` a one-to-one bare `Allele` with no `Variant` at all, on the grounds that a fusion has no
coordinate. That buys the classification chain and nothing else, and it costs the two things that
matter most for somatic analysis:

- Gene-list filtering and compound-het two-hit detection both join `VariantGeneOverlap`, which is a
  `Variant` FK — `snpdb/variant_queries.py:31`, `analysis/models/nodes/sources/trio_node.py:221-232`,
  `quad_node.py:190-198`. No `Variant` means fusions are structurally invisible to gene lists and
  comp-het, not merely absent from a grid.
- Sample-level presence *is* `CohortGenotype.variant` (`snpdb/models/models_cohort.py:733`). Analysis
  is sample-scoped, so with no `Variant` a fusion cannot appear in any analysis node.

**The file has breakpoints.** All 33 rows of the test file carry `Gene A Breakpoint` and
`Gene B Breakpoint` (`chr8:128806980` form), 12 of them inter-chromosomal, with
`Fusion Directionality Known=True` and per-side sense on every row. The "gene-level fusions may have no
base-level coordinates" question in the issue body is answered for TSO 500: they do. What is deferred is
the *representation* of breakends, not the acquisition of the data — capture the values from day one.

**Coordinates are per-observation, so identity stays gene-level.** In the one test file:

```
ENTPD3-RPL14  DRAGEN  chr3:40442308 → chr3:40503522
ENTPD3-RPL14  DRAGEN  chr3:40446552 → chr3:40503522
ENTPD3-RPL14  DRAGEN  chr3:40447887 → chr3:40503522
```

One caller, one sample, one gene pair, three 5′ breakpoints. Coordinate-keyed identity would make that
one biological event into three fusions — and this is *within* a caller, before the DRAGEN-vs-SpliceGirl
comparison. So the gene pair is the identity and the breakpoints are observation data.

## Storage

### The gene-level contig

Gene-level events go on **one contig of their own, shared by every genome build**, rather than on the 5′
gene's real chromosome at a synthetic position.

The alternative — anchoring on the real chromosome with `position = 1_000_000_000 + HGNC ID` — fails the
contig narrowing that analysis applies. `AnalysisNode.get_queryset` adds
`Q(locus__contig__in=self.get_contigs())` (`analysis/models/nodes/analysis_node.py:651-652`) and the
gene-derived nodes build that set from transcripts (`gene_list_node.py:91-99`, `moi_node.py:173-176`,
`phenotype_node.py:137-140`). A fusion's `VariantGeneOverlap` rows name *both* partners, so for any
inter-chromosomal fusion the 3′ gene's contig differs from the anchor's, and a gene list on that partner
silently drops it — EML4-ALK and EGFR-SEPTIN14 pass, CD74-ROS1, ETV6-NTRK3 and BCR-ABL1 fail. Half the
canonical set working is worse than none of it working.

Pair contigs (`chr9_chr22`) were considered and rejected: they put chromosome assignment back inside
`Variant` identity — the exact instability that keying on HGNC ID avoids — the pair space is not 24×24
once scaffold and clone-based partners are counted (`RP11-458D21.5`, `AC016683.6` are both in the test
file), they need a lookup model to be queryable, and single-gene CN events have no pair.

One contig, then:

- `name` obviously non-genomic (`GENE_LEVEL`), `refseq_accession` a VG-namespaced string so it cannot
  collide with a real accession.
- A new `SequenceRole` value, so the `role=ASSEMBLED_MOLECULE` filter in `standard_contigs`
  (`snpdb/models/models_genome.py:147-148`) keeps it out of the chrom mappings that VCF import uses.
- `length` set well above any future HGNC ID, since `Variant.validate` bounds position by it
  (`snpdb/models/models_variant.py:650-659`).
- A `GenomeBuildContig` row for every build. `Contig` deliberately does not link to a build so builds can
  share them (`models_genome.py:379-382`), and liftover already exploits that: `SAME_CONTIG` /
  `_liftover_using_existing_contig` is the "Mito, 37 and 38 contigs are the same so we can re-use a
  variant" path (`snpdb/liftover.py:239-252`, `AlleleConversionTool.SAME_CONTIG` at
  `snpdb/models/models_enums.py:158`). So one `Variant` serves every build and liftover needs no new
  conversion tool and no coordinate translation.

Consequence worth stating plainly: the fusion `Variant` is build-independent, and "which build the reads
were aligned to" lives on the VCF/Sample, exactly as it does for every other variant. The alternative —
one gene-level contig per build — is available if that ever needs asserting at variant level, at the cost
of a custom liftover tool.

`GenomeFastaContig` will have no row for this contig, so anything asking for its fasta name raises
`ContigNotInFastaError` (`models_genome.py:490`). That is a per-lookup failure, and gene-level variants
are excluded from every VCF-writing path below — but confirm it before the migration lands.

### Locus and Variant

- `Locus.contig` = the gene-level contig, `Locus.position` = the anchor gene's HGNC ID,
  `Locus.ref` = `N`.
- `Variant.end` = `position`, `Variant.svlen` = `0`.

`svlen = 0` rather than null is deliberate: `unique_together = ("locus", "alt", "svlen")`
(`models_variant.py:592-593`) is a Postgres unique constraint, and Postgres treats nulls as distinct, so
a null `svlen` gives no database-level uniqueness. The cost is that `Variant.get_symbolic_q()` (svlen not
null, `:637-638`) then calls gene-level events symbolic, which routes them at the SV annotation pipeline
— they need excluding there regardless, since VEP cannot parse the alt, so this adds no new work.

### The identifier space — `FusionGeneId`

Anchors and partners are named by a `FusionGeneId` row rather than straight off `HGNC`:

```python
class FusionGeneId(models.Model):
    CUSTOM_ID_START = 1_000_000
    symbol_str = TextField(unique=True, db_collation='case_insensitive')
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)
```

`pk` **is** the HGNC ID where the gene has one, so an ordinary fusion's identity is the same number on
every deployment. Symbols HGNC does not carry get a pk above `CUSTOM_ID_START`.

That last part is why the table exists. Clone-based identifiers are routine fusion partners — the test
file has `RP11-458D21.5;NOTCH2NL` and `PPARG/AC016683.6` — and keying identity on HGNC alone leaves
them unrepresentable, forcing a park-or-guess decision on rows the caller reported perfectly clearly.
With a local id every call becomes a `Variant`, and `<FUSION:UNKNOWN>` goes back to meaning only what
it says: the caller named one gene and left the other unspecified.

**Local ids do not travel.** Anything leaving the system — display, export, a classification sent to
another instance — uses `symbol_str`, and the receiver resolves it through its own table.
`GeneFusion.canonical_str` is that string form. This is stated on the model, because the numbers look
portable and half of them are.

### Identity in the alt

Biological identity is encoded in the alt string, which hashes to its own `Sequence` row and therefore
its own `Variant`:

- `<FUSION:HGNC:nnn>` — directional. Anchor is the 5′ partner; the alt carries the 3′ partner's ID.
- `<FUSION_UNORDERED:HGNC:nnn>` — both partners known, direction not asserted. Anchor is the gene with
  the smaller ID; the alt carries the larger.
- `<FUSION:GENE:nnn>` / `<FUSION_UNORDERED:GENE:nnn>` — same, for a partner with no HGNC entry. The
  namespace says whether the number travels, so nobody has to know that 1,000,000 is the dividing line.
- `<FUSION:UNKNOWN>` — the caller named one gene (the anchor) and no partner at all.
- `<FUSION_UNORDERED:UNKNOWN>` — one gene named, and it is not the 5′ side. Anchoring a 3′ partner as
  if it were 5′ would assert a direction the caller did not.

Reciprocal fusions are distinct variants, because the 5′ promoter drives a different protein and the
clinical significance differs (BCR-ABL1 vs ABL1-BCR is the canonical case). The test file settles this
rather than theory: `ENTPD3→RPL14` appears three times and `RPL14→ENTPD3` once, same caller, same
sample. Ordered and unordered are distinct, because an unordered report does not assert direction.

HGNC-backed identity is preferred wherever available: HGNC bridges RefSeq and Ensembl, so it stays
consortium-agnostic, and HGNC IDs do not move, which symbols demonstrably do — the test file's `SEPT14`
is HGNC's `SEPTIN14`, and `cnv.vcf`'s `SEGID` says `MYCL1` where the rest of the pipeline says `MYCL`.
Delimiter normalisation happens at the import layer and never reaches identity: `;` and `/` separate
names *within* one side's cell (a hyphen cannot — `RP11-458D21.5` contains one), while `::`, `~`, `--`
and a lone `-` separate the two sides of a fusion written as one string.

### `GeneFusion` companion model

One-to-one with `Variant`, denormalising what grids and filters need to join and sort on:

```python
class GeneFusion(models.Model):
    variant = models.OneToOneField(Variant, on_delete=CASCADE)
    anchor = models.ForeignKey(FusionGeneId, related_name='fusions_as_anchor', on_delete=PROTECT)
    partner = models.ForeignKey(FusionGeneId, null=True, related_name='fusions_as_partner', on_delete=PROTECT)
    is_ordered = models.BooleanField(default=False)
```

`partner` is null only for the UNKNOWN forms. `is_ordered` is derivable from the alt prefix but
denormalised for filtering. `clean()` enforces that the alt prefix, the alt's gene ID and the two FKs
agree.

The `Allele` arrives the ordinary way, via `VariantAllele`, rather than being owned by `GeneFusion` —
one allele per fusion, build-independent, with the single shared-contig variant beneath it.

**No `source_data` blob on this model.** The earlier design put breakpoints, exon designations and the
original fusion string there, which the three-row `ENTPD3-RPL14` case rules out: those values belong to
the call, not to the fusion. They go in `CohortGenotype.info` with the rest of the per-observation data.

## Containment

Every place that has to know a gene-level variant has no coordinate, and what it becomes:

| Where | Today | Change |
|---|---|---|
| `annotation/annotation_version_querysets.py:26-39` | `pipeline_type_variant_q()` — documented single source of truth for which variants a pipeline claims | Gene-level predicate registers here; claims neither STANDARD nor STRUCTURAL_VARIANT |
| `annotation/models/models.py:1829-1833` | `VARIANT_ANNOTATION_Q` — "filters to describe variants that can be annotated" | Exclude gene-level |
| `snpdb/models/models_variant.py:783-802` | `clingen_allele_skip_reason()`, `can_have_c_hgvs`, `can_have_annotation` | Gene-level returns a skip reason; no ClinGen call, no c.HGVS attempt |
| `snpdb/models/models_variant.py:398-419` | `as_external_explicit()` raises `Unknown symbolic alt` at `:415` | Guard beside `can_be_made_explicit` |
| `analysis/models/nodes/analysis_node.py:494-514, 651-652` | Contig narrowing intersects gene-derived contig sets | `get_contigs()` unions in the gene-level contig wherever it narrows — one place, inherited by every gene-derived node |

A named `Variant.get_gene_level_q()` (and its instance-level twin) is the single predicate all of these
call, so the special case is one concept rather than scattered checks, and review can grep for it.

The contig narrowing is gated on `Analysis.node_queryset_filter_contigs`, default `False`
(`analysis/models/models_analysis.py:63`), so it is an opt-in optimisation rather than the default path —
but it needs to be correct before any deployment turns it on.

## Annotation — a third pipeline type

Gene-level variants never reach VEP. Their annotation is computed locally from the gene identity, as a
new `VariantAnnotationPipelineType.GENE_LEVEL` (`annotation/models/models_enums.py:77-80`).

A pipeline type rather than a one-off at import, because:

- `VariantGeneOverlap` is keyed `(version, annotation_run, variant, gene)`
  (`annotation/models/models.py:2228-2236`) and symbol-to-gene resolution is per `GeneAnnotationRelease`,
  so the rows must be rebuilt for each new `VariantAnnotationVersion` exactly as VEP output is. Written
  once at import, fusions would silently drop out of gene lists at the next annotation version.
- `VariantGeneOverlap.annotation_run` is non-null, so the rows need a run to belong to.
- `AnnotationRun.get_for_variant` picks its pipeline type from `variant.is_symbolic`
  (`annotation/models/models.py:1145-1156`) and wants the third branch.

What a run writes per fusion:

- `VariantGeneOverlap` rows for **both** partners. This is what makes a gene list on either side find the
  fusion, and what lets a fusion count as one of the two hits in comp-het.
- A `VariantAnnotation` row with `gene_symbol` and `overlapping_symbols` populated, so grid columns and
  hover cards have something to show. `annotation/management/commands/fix_annotation_sv_overlaps.py` and
  `SVGeneOverlapResolver` are the shape to copy — the same backfill for long SVs VEP skipped as
  `TOO_LONG` — minus the coordinate resolver, since the genes are known outright.
- A new `VEPSkippedReason` value (`annotation/models/models_enums.py:150-154`) marking the row as locally
  computed, so an empty consequence field is not mistaken for "VEP found nothing".

Resolution path is `FusionGeneId` → symbol → release genes, not HGNC → gene. `HGNC` carries
`ensembl_gene_id` but no Entrez ID (`genes/models.py:77-97`) while `Gene.identifier` is the Entrez ID
for RefSeq releases, so RefSeq resolves through the symbol and the per-release symbol cache
(`ReleaseGeneSymbol` / `ReleaseGeneSymbolGene`) — the same resolver Phase 6 of the overall plan needs
for `cnv.vcf`'s `SEGID` symbols. Going through `FusionGeneId` rather than `HGNC` also means a
locally-identified partner resolves by exactly the same path as an HGNC-backed one.

Identity keys on the id for stability and resolves through the symbol at annotation time. That is the
right split: identity is fixed forever, resolution is versioned and re-runnable, so a later HGNC import
that fixes `SEPT14` → `SEPTIN14` improves the mapping without touching a single `Variant`.

**How it runs.** `_handle_range_lock` already creates one `AnnotationRun` per member of
`VariantAnnotationPipelineType`, and the dispatcher already launches a CREATED run through
`annotate_variants`, so the scheduler needs no change at all. `annotate_variants` branches on
`pipeline_type` to call `annotate_gene_level_run` instead of `dump_and_annotate_variants` — the same
shape as STANDARD and STRUCTURAL_VARIANT sharing that task and branching inside it. Everything around
the branch (the task_id lock, the blocker check, the lease heartbeat, error reporting, releasing the
lease and kicking the dispatcher) is identical for all three types, so it stays in one place.

The one thing gene-level does differently is finish there rather than at ANNOTATION_COMPLETED: with no
VEP subprocess to throttle and only a handful of rows to write, there is nothing for the second lane
(`db_workers` → `import_annotation_run`) to do, so `annotate_gene_level_run` writes the rows, walks the
status fields to FINISHED and sends `annotation_run_complete_signal` itself. Gene-level variants are
rare, so `count_annotation_run` finishes almost every GENE_LEVEL run as empty without dispatching it.

### Why not feed VEP a BND pair

The breakpoints exist, so BND records are writable. They still are not worth writing:

- Whether VEP parses BND ALT syntax at all is unverified — worth a 20-minute experiment before anyone
  designs around either answer.
- Parsed or not, VEP's output is per-position feature overlap: which transcript and exon a breakend
  lands in. It does not do frame analysis, retained domains, or fusion recognition, which is the
  annotation that would actually be worth having (AGFusion/Arriba territory).
- The part VEP would give us is available twice over already: computable locally from the release's exon
  structures, and reported by the caller as `Gene A Location` / `Gene B Location`
  (`IntactExon`/`BrokenExon`/`Intronic`/`Intergenic`).
- BND is two records, so two `Variant`s. `Allele.variant_for_build_optional` resolves through `.first()`
  (`snpdb/models/models_variant.py:126-134`), so a second `VariantAllele` for the same allele and build
  is not rejected — it is silently ignored, and which variant you get depends on ordering. One variant
  per allele per build is effectively an invariant of the allele layer.

## Ingestion — the `AllFusions.csv` parser

An import task factory claiming the file, recognisable from the `# Source = FusionProcessor` first line
and the `Caller,Gene A,Gene B,…` header row, with a processing ability above `GeneListImportTaskFactory`'s
default 1 (`upload/import_task_factories/import_task_factory.py:35-40`), which would otherwise win and
import the fusions as a gene list.

**Through the normal VCF import pipeline, minus the bcftools stages.** The rows are written out as a
VCF of gene-level variants and inserted by the standard path, so the insertion queue, `VariantPKLookup`
and per-pipeline-type max-variant tracking all apply unchanged - `UploadedVCFPipelineMaxVariant` is
keyed on `VariantAnnotationPipelineType` and populated via `pipeline_type_variant_q`, so GENE_LEVEL
tracking (and therefore the annotation gate) needs no new code.

What is skipped is only what needs a reference base: `bcftools norm --check-ref=s` reads it from the
fasta, and `vcf_clean_and_filter` maps chroms through `standard_contigs_only` (which excludes the
gene-level contig by role) before replacing the header with one built the same way. Nothing is lost -
the VCF is written by us rather than supplied by a lab, so it already has our contig names, no
undeclared FILTERs, one alt per record and sorted order. `preprocess_gene_level_vcf` therefore goes
straight to the split, and everything downstream of that is the shared code path.

`END = POS` on each record is what gives `svlen` 0 through
`vcf_get_ref_alt_svlen_and_modification`, which requires SVLEN or END for any symbolic alt.

**The file creates its own `VCF` and `Sample`.** A `Sample` belongs to exactly one `VCF`
(`snpdb/models/models_vcf.py:334`) and `CohortGenotype` hangs off a per-VCF collection, so fusion rows
cannot be appended to the RNA arm's existing `SpliceVariants` sample. `VCF` has no file FK
(`snpdb/models/models_vcf.py:68-102`), so creating one for a CSV is cheap. Phase 4 of the overall plan is
what ties that sample to the extraction; Phase 7's grouping node is what shows both arms together.

**Per-observation data goes in `CohortGenotype.info`** (JSONField, `models_cohort.py:751`) — caller, both
breakpoints, score, filter string, split/pair read counts, `Gene A/B Location`, `Gene A/B Sense`, and the
fusion string as written.

**Record the collapse.** Several rows can share one gene pair (the `ENTPD3-RPL14` case), and they become
one `Variant` and therefore one `CohortGenotype` row. That is the same merge Phase 2 of the overall plan
built `ModifiedImportedVariantOperation.SHARED_LOCUS` for (`upload/models/models_enums.py:56`), with
`operation_detail` carrying what was summed (`upload/models/models.py:729-730`) so per-row values stay
reconstructable. Use the same discipline: keep every row's data in the info list, and record that a merge
happened. The loader creates one `UploadStep` for the import, which is what `ModifiedImportedVariants`
hangs off.

**The genome build is declared, not detected.** The file's breakpoints are `chrN:pos` with no header to
detect a build from, so a multi-build deployment needs `genome_build` upload metadata — the same key
`_DragenExonCNV.vcf` already requires. The fusion `Variant` itself is build-independent, so this is only
for the VCF/Sample bookkeeping and the `VariantAllele`.

**Ingest unfiltered.** 1 of 149 rows in the real run passed the caller's own filter, and the filter
strings are long semicolon-joined lists. Same posture as the CNV and small-variant VCFs: take everything,
apply our own thresholds. Carry `Caller` per row — two callers write into one file and the header warns
their scores and filters are not comparable.

**Row shapes to survive**, all present in the test data:

- Multi-gene partners are routine, not edge cases — `ROS1;GOPC`, `FIP1L1;PDGFRA`,
  `RP11-458D21.5;NOTCH2NL` — with the per-gene annotation columns correspondingly semicolon-delimited.
- Mixed delimiters — `PPARG/AC016683.6` uses a slash alongside the `;` and `-` forms.
- `SEPT14` where HGNC says `SEPTIN14`, and spreadsheets say a date.
- Direction reported rather than inferred — `Fusion Directionality Known` plus `Gene A/B Sense` populate
  `is_ordered` from data.

## Classification linkage

`ImportedAlleleInfo` gains a `gene_fusion` FK and a short-circuit that skips HGVS resolution: with the
fusion known, the allele is `gene_fusion.variant`'s allele and the status is matched immediately. A lab
submitting `BCR::ABL1` as a classification target resolves through the same canonicalisation the parser
uses.

Nothing else in the classification chain changes. `ClinicalContext`, `DiscordanceReport` and ClinVar
export all key on `Allele`, which fusions now have in the ordinary way.

## Gene-level CNV — designed, not built

TSO 500 v2.6.2's `cnv.vcf` emits `<DUP>`/`<DEL>` with real coordinates and `FORMAT/SM`, and it already
imports. Those stay coordinate-based structural variants. `<AMP:HGNC:nnn>` / `<LOSS:HGNC:nnn>` on the
gene-level contig is only for callers that report a gene-level event with no coordinates at all — the
`CombinedVariantOutput` "JAK2 amplification (5 copies)" style. Same anchor, same contig, same annotation
run; copy number stays observation-level in `CohortGenotype.info` rather than in the alt, because
different labs use different amplification thresholds and per-count identity would stop two labs ever
agreeing on "JAK2 amplification".

Design it now, build it when a file needs it. It is purely additive — one enum value and one alt prefix.

## Order of work

1. The gene-level contig: `SequenceRole` value, `Contig` + `GenomeBuildContig` migration, and the
   `Variant.get_gene_level_q()` predicate.
2. Containment — the rows of the table above, each with a test.
3. `FusionGeneId`, the `GeneFusion` model, canonicalisation (delimiters, alias resolution,
   ordered/unordered/unknown), and `get_or_create` from a gene pair.
4. The `GENE_LEVEL` annotation pipeline type and its run.
5. The `AllFusions.csv` factory and parser: VCF + Sample creation, direct variant insert, `CohortGenotype`
   with per-row info and merge recording.
6. `ImportedAlleleInfo.gene_fusion` short-circuit, and the search receiver.

Steps 1-2 are the foundation everything else sits on. Step 4 is independently testable against
hand-created fusion variants, so it can run in parallel with step 5.

Note on `VARIANT_ANNOTATION_Q`: the containment table above lists it as excluding gene-level, and it
does not. It means "variants that get a `VariantAnnotation` row", which gene-level variants do — just
written by the GENE_LEVEL run rather than VEP — and `get_variants_qs_for_annotation` builds every
pipeline's queryset from it, so excluding there would leave the GENE_LEVEL pipeline with nothing to
annotate. The exclusion belongs in `pipeline_type_variant_q` alone, which is the documented single
source of truth and subtracts gene-level from both VEP types.

## Search

A fusion has to be findable by the name people call it: `BCR::ABL1` and `CD74-ROS1` both resolve
through a `Variant` search receiver, either ordering, with aliases applied (`EGFR::SEPT14` finds
`EGFR-SEPTIN14`). Search is **lookup only** — it runs on whatever a user types, so it must never mint a
`GeneFusion` or a `FusionGeneId`, unlike the loader and the classification short-circuit, which both
create on demand.

Two things this leaves for later: finding a fusion by a single partner (`ROS1` returning every fusion
it takes part in) belongs with Phase 6's grid work, and searching a `GENE:`-namespaced local id is
deliberately not offered — the number means nothing to a person and nothing on another deployment.

## Decisions still open

| Decision | Why it matters |
|---|---|
| Does VEP parse BND at all? | Only affects whether a future breakend phase can reuse VEP; nothing in phase 1 waits on it |
| Single grid or multiple grids for non-variants | Overall plan Phase 6, best answered against real fusion rows |

Settled during implementation: **a multi-gene side resolves to its first HGNC-backed member**, falling
back to its first known gene symbol, then to the first name as written. A side that resolves to nothing
known still gets an identity through `FusionGeneId`'s local id space, so there is no park-or-skip case
and no row is lost. The cell is kept exactly as the caller wrote it in `CohortGenotype.info`.

## Done when

All 33 rows of `ExampleSample_RNA_2600000001B_AllFusions.csv` import; `EGFR-SEPTIN14` resolves despite
the file saying `SEPT14`; the same gene pair from both callers — and the three `ENTPD3-RPL14` rows from
one caller — resolve to one `GeneFusion` with every row's breakpoints preserved in `CohortGenotype.info`;
a gene list containing `ROS1` finds `CD74-ROS1` with contig narrowing turned on; searching `BCR::ABL1`
finds the fusion; and the fusion sample's variants carry gene symbols on the grid.

## Deferred, deliberately

- **Fusion equivalence and discordance** — the research-level part of #1506.
- **Breakend representation** — a per-observation breakend model, and BND VCF export. The values are
  captured from day one, so this is a read-side change when a consumer appears.
- **Frame/domain fusion annotation** — not something VEP does; a separate tool integration if wanted.
- **Ensembl-only or RefSeq-only partners** — `<FUSION:ENSG:nnn>` is additive if evidence accumulates.
- **Gene-level CNV loading** — designed above, built when a coordinate-free caller arrives.
