# TSO 500 — order of work

Sequencing plan for [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431)
and the issues under it. Covers what to build in what order and why, not how to build each piece —
each issue carries its own design.

The file-level ingestion analysis is in [`tso500_ingestion_plan.md`](tso500_ingestion_plan.md); the
test data and its gotchas are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).
Three phases still ahead have their own design docs:
[`fusion_variants_issue_1506_plan.md`](fusion_variants_issue_1506_plan.md) for gene fusions,
[`tso500_phase_7_plan.md`](tso500_phase_7_plan.md) for the grouping node, and
[`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md) for Phase 8's
reporting. Phases 2, 3 and 4 had docs of their own, deleted once they landed — what outlived them is in
the Done entries below, and the PRs that reference them are #1712, #1715 and #1716.

---

## Done

- **Phase 0 — pipeline fixes** (PR #1709, branch `issue_1506_tso500_pipeline_fixes`).
  `vcf_header_filter_ids()` in `library/genomics/vcf_utils.py` reads declared FILTER IDs off the raw
  header lines, so `_DragenExonCNV.vcf` no longer dies at the first pipe stage on PyVCF's FILTER regex.
  A regex rather than cyvcf2, which needs a real file descriptor and so does not suit a stdin filter.
  Also stopped `PatientRecordsImportTaskFactory.get_processing_ability` raising out of file-type
  detection on a CSV pandas cannot read, which was breaking detection for every other CSV type.

- **Phase 1 — `Specimen → Extraction → Sample` (#1704)** (PR #1705, branch
  `issue_1704_specimen_extraction_sample`). The join key everything below groups on:

  ```
  Patient └── Specimen └── Extraction └── [SequencingSample] └── Sample
  ```

  `nucleic_acid_source` moved off `Specimen` onto the new `Extraction`, so the DNA and RNA arms of one
  block share a specimen instead of becoming two. `Specimen` swapped its `TextField` primary key for a
  surrogate one with `reference_id` unique per patient. Both levels extend `ExternallyManagedModel` and
  carry a local `reference_id` beside the nullable `external_pk`, since not every deployment has a system
  managing its records. `Sample.extraction` replaced `Sample.specimen`; `SequencingSample.extraction`
  added as optional enrichment. R7's consumers updated, the patient CSV still round-trips through one
  extraction per specimen, and the patient page gained an Extractions tab. Everything in `patients` is now
  a `TimeStampedModel`, and `Extraction.extraction_date` tells re-extractions apart.

- **Phase 2 — per-file loader refinements** (PR #1712, branch
  `issue_1711_tso500_phase2_loader_refinements`, issue #1711). All four items:

  `FileUpload.metadata` is a JSON blob carrying the facts a file doesn't reliably hold itself, with each
  `ImportTaskFactory` declaring the keys it accepts and validation running at upload time — so a mistyped
  key is a 400 while the client is still connected rather than a `REQUIRES_USER_INPUT` several stages
  later. A declared `genome_build` that contradicts header detection fails the import rather than picking
  a winner. `VCFSourceSettings` grew a `sample_field_overrides` JSONField (one column, so "clear this
  field" is expressible), and the `^SpliceGirl` row remaps `AD`/`DP` so VAF derives as
  `ALTDEDUP/(ALTDEDUP+REFDEDUP)`. FILTER values the source header never declared are moved into
  `INFO/VG_UNDECLARED_FILTERS` by `vcf_clean_and_filter` and restored by the genotype processor at
  insert — they can't stay in the FILTER column, because `bcftools norm` *dies* on an undeclared FILTER
  rather than warning, and dropping them loses real calls. Copy-neutral records
  (`ALT=.` + `END` + no `SVTYPE`) are skipped and counted;
  `SVTYPE` marks a caller describing an event rather than a segmentation interval, so `_DragenExonCNV.vcf`'s
  per-gene `Undetermined` call survives. Where separate VCF rows share a locus and have their depths
  summed, a `ModifiedImportedVariant` with the new `SHARED_LOCUS` operation records the summed depths so
  the per-record VAF stays reconstructable.

- **Phase 3 — tissue status, and pages for Specimen / Extraction** (PR #1715, branch
  `issue_1706_tso500_phase3_tissue_status_and_pages`, issues #1706 and private#2447). Both Phase 1
  leftovers, done before Phase 4 puts `Specimen` behind a public API:

  `Specimen.mutation_type` became `Specimen.tissue_status` — Reference / Affected / Unknown, not null,
  defaulting to Unknown. The old field was Germline/Somatic with `default=GERMLINE`, so every tumour
  block nobody had touched was stamping its classifications "Germline" through autopopulate and filing
  them alongside other labs' real germline records. It sits on `Specimen` rather than on `Tissue`
  because the same tissue plays different roles — blood is the reference in a solid-tumour workup and
  *is* the tumour in leukaemia.

  Three levels answer three different questions, and the phase settled which field owns each — worth
  quoting rather than rederiving, since #1707 serializes the first and #1559 and Phase 8 read all three:

  | Level | Question | Where | Values |
  |---|---|---|---|
  | Specimen | What is this material, and what role does it play in the test? | `Specimen.tissue_status` | Reference / Affected / Unknown |
  | Call set | What's in this VCF? | `Sample.variants_type` | Unknown / Germline / Mixed / Somatic only |
  | Variant | What is *this* variant's origin? | `allele_origin` → `AlleleOriginBucket` | germline / somatic / other |

  private#2447 says `VCF.variants_type`; it is actually on **`Sample`**
  (`snpdb/models/models_vcf.py:343`), alongside `Sample.extraction` and `Sample.is_somatic`. That is why
  the derivation sits in `get_evidence_fields_for_sample_and_patient` rather than the extraction helper,
  which cannot see the call set.

  Autopopulate became a derivation across two of those levels instead of one field read: an
  `allele_origin` is asserted only where the specimen and `Sample.variants_type` agree — reference +
  germline is `germline`, anything + somatic-only is `somatic`, everything else stays unset and falls
  back to `settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET`, because a mixed tumour sample's origin is
  genuinely per-variant and cannot be known at accession. It emits the evidence key's own option string
  rather than a display label, which is what `bucket_for_allele_origin` matches on.

  `tissue_status` is per-specimen only — both arms of one block share the material, so it did not go on
  `Extraction` the way `nucleic_acid_source` did. `Tissue` stays a lookup table with no page of its own.
  Matched normal stays derivable — "reference specimen, same patient" — rather than becoming an FK.

  `patients/0016` maps `S` → Affected and `G` → Unknown, `G` being indistinguishable from an untouched
  record, and registers a `ManualOperation` only where a `PatientRecord` actually carries a `G`. The
  patient CSV column is a hard rename to `'Specimen Tissue status (Reference/Affected/Unknown)'` — an
  old spreadsheet fails the import naming the missing column, since Germline → Unknown is a loss of
  meaning rather than a translation.

  Both models now delegate permissions to their patient the way `Trio` delegates to its cohort, which
  gave them `filter_for_user` and let the autocompletes drop their hand-rolled patient filters. Detail
  pages at `view_specimen` / `view_extraction`, editable except `external_pk`; samples on both pages
  filter through `Sample.filter_for_user` separately, since a sample carries its VCF's permissions
  rather than its patient's. `PreviewModelMixin` on both, search on `HAS_3_ANY` (a TSO 500 reference is
  all digits), `ExternalPK` search reaching both, and links each way along
  `Patient → Specimen → Extraction → Sample`.

- **Phase 4 — how identifiers cross the lab boundary** (PR #1716, branch
  `issue_1707_tso500_phase4_patient_api_and_extraction_linking`, issues #1707 and #1559). #431's
  decision governs it — **upload the files separately and join later** — so the whole phase is
  assignment on top of the per-file upload that already existed, in four steps: resolve, create, name,
  reconcile, then measures on top.

  `patients/external_references.py` is the one place a client's identifier becomes a row. A bare string
  is the local reference (`reference_id` / `patient_code`); reaching an `ExternalPK` takes `code` **and**
  `external_type` together, since `ExternalPK` is unique on `(code, external_type, external_manager)` and
  a code alone names nothing — rejected on shape while the client is still connected. Resolution has
  three outcomes rather than two: nothing found is **Pending**, because an extraction legitimately
  arrives after the VCF that names it, and more than one found is **Needs attention** — `reference_id` is
  unique per parent rather than globally, so that is ambiguity rather than a precedence rule.

  #1707's API (`api/v1/patient`, `specimen`, `extraction`, `specimen_measure` + `bulk_create`) upserts
  on the identifiers sent, so a client re-posting a run gets the same rows back. A specimen naming an
  unknown patient is a 400 — it has nowhere to live — unlike a VCF naming an extraction, which has a row
  to park on. An unknown `external_manager` is a 400 naming the configured ones unless the poster is a
  superuser (`PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY`), since creating one silently would decide
  `can_modify` off a misspelling. `/patients/api/` joins `PUBLIC_PATHS` because
  `GlobalLoginRequiredMiddleware` runs before DRF and would redirect API calls; `IsAuthenticated` still
  applies.

  Two routes let a VCF name its extraction: `POST api/v1/sequencing_sample/link_extraction`, one call per
  sequencing sample reaching all of that arm's VCFs through `link_samples_and_vcfs_to_sequencing` and
  carried across when a sheet is re-sent; and `extraction` / `sample_extractions` upload metadata keys,
  which need no seqauto records at all. Two keys rather than one so each carries a single type and VCF
  sample names share a namespace with nothing. Where the two routes name different extractions the import
  fails, the same rule a declared `genome_build` the header contradicts gets.

  `ExtractionMatchMixin` parks the claim beside the FK on both `Sample` and `SequencingSample`, so an
  unresolvable reference is parked rather than rejected. `reconcile_pending_extractions` re-resolves
  hourly and again whenever an extraction is created, promoting stale Pending to Needs attention after
  `PATIENT_EXTRACTION_MATCH_PENDING_DAYS = 3`. `extraction_match_date` moves only when the claim itself
  changes — stamping it every pass renewed the clock, so a Pending row could never age out. The task also
  carries an extraction down to `Sample`s whose link call arrived *after* their VCF imported, which
  nothing else can, `link_samples_and_vcfs_to_sequencing` running once at import.
  `PATIENT_EXTRACTION_SAMPLE_NAME_REGEX` derives a reference from the VCF sample name for deployments
  with nothing upstream to quote one — off by default, consulted only where nothing was posted, marked
  `derived` on the parked row, and it creates no `Specimen` or `Extraction`.

  #1559's `SpecimenMeasure` is vendor-neutral and client-posted rather than parsed out of pipeline
  output, keeping both the score and the call plus the raw payload, one current value per
  `(specimen, measure_type)`, shown on `view_specimen`. Parked references surface on `view_sample`, in
  the `Sample` / `SequencingSample` admin filters, and through a health check that stays silent while
  everything matches.

## Where things stand

The five test files from `ed5e15a33` — one de-identified run, one specimen, DNA and RNA arms:

| File | State |
|---|---|
| `hard-filtered.vcf` | loads, 93 records |
| `cnv.vcf` | loads, 16 records — 9 copy-neutral rows skipped and counted; `SM` surfacing is Phase 6 |
| `SpliceVariants.vcf` | loads, 17 records — VAF derived from `ALTDEDUP`/`REFDEDUP`, `LowUniqueAlignments` preserved |
| `_DragenExonCNV.vcf` | loads, 2 records, against a `genome_build` declared at upload |
| `AllFusions.csv` | no longer breaks file-type detection; needs a parser and somewhere to put the rows (Phase 5) |

## Dependency map

```
  Phase 0  pipeline fixes ✔ ──► Phase 2  loader refinements ✔
           chase TAU files ─────────────────────────────────────────┐
                                                                    │
  Phase 1  #1704 Specimen→Extraction ✔                              │
             │                                                      │
      ┌──────┴───────────────────────────┐                          │
      │                                  │                          │
  P3  private#2447 tissue status ✔       │                          │
      #1706 pages, search, preview ✔     │                          │
      │              │                   │                          │
  P4  seqauto link ✔ ┤                   │                          │
      #1707 patient API ✔ ► matching ✔ ► #1559 measures ✔ ◄─────────┘
      │              │                   │
      │              │              P7  private#223 grouping node
      │              │
      │         P8  sapath#246 ──► #444 multi-variant reporting
      │
  P5  #1506 GeneFusion ──► AllFusions parser   (independent of #1704 — can run in parallel)
      │
  P6  #1558 grid display  (+ #1706's grid half, private#2837 — needs P3's URLs)
```

Phases are listed in the order to do them, but only the arrows are real constraints. #1506 is the one
remaining critical-path item — the fusion parser has nowhere to put rows without it. Everything else is
ordered by value.

---

## Still open from Phase 0 — chase TAU for the missing files

Long lead time, so worth pushing now:

- `Logs_Intermediates/Gis/<sample>/<sample>.abcn_annotated.vcf` + `<sample>.abcn_genes.tsv` — per-gene
  absolute and minor copy number, the latter being where gene-level LOH comes from. The only genuinely
  CVO-only data, and a VCF, so it is in ingestion scope. Illumina names the path but publishes no
  INFO/FORMAT spec, so it cannot be mocked up.
- `MetricsOutput.tsv` / the run-level TMB summary — where the numbers Phase 4's `SpecimenMeasure` API
  takes are transcribed from, so it is now the client's need rather than VG's.
- A run with a real BRCA1/2 large rearrangement. Both `_DragenExonCNV.vcf` records in the test data are
  constructed, and no published TSO 500 output has a populated one either. Until one arrives the field
  spelling in that file stays provisional.

## Still open from Phase 1

The patient CSV columns are all `SPECIMEN_*`, so a CSV can name one extraction per specimen but not the
DNA and RNA arms of one block. Worth revisiting when a consumer needs it.

## Still open from Phase 2

`AllFusions.csv` still has no parser — that is Phase 5. The manual checklist on #1711 wants running
against a real GRCh37 deployment: PR #1712 measured the record counts and VAFs through
`write_cleaned_vcf_header` → `vcf_clean_and_filter` → `bcftools norm`, not through a full import with
annotation.

The `source` strings a client sends (`DRAGEN TSO500 SmallVariant`, `DRAGEN TSO500 CNV`) are documented
in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md) but nothing is keyed
on them yet — the SpliceGirl mapping comes off the header and the copy-neutral skip is a general rule.
They become part of VG's configuration contract the moment a `VCFSourceSettings.source_regex` matches
one, so they want to stay stable from the first client.

`SEGID=MYCL1` resolving to `MYCL` is unconfirmed against a real database. NCBI carries the alias, and
`GeneSymbolMatcher.get_gene_symbol_id_and_alias_id` (`genes/gene_matching.py:45`) is the resolver, but
Phase 6's gene-symbol item rests on that assumption holding — and Phase 5's fusion parser wants the same
shape for `SEPT14` → `SEPTIN14`, so one check covers both.

**Clients send a build's own name (`GRCh37`), not an alias.** `GenomeBuild.get_name_or_alias("hg19")`
raises `MultipleObjectsReturned` rather than `DoesNotExist`, so a declared build that will not resolve
is a 400 rather than a guess. `hg19` happens to resolve here because that build is disabled, but
`GenomeBuild.enabled` is per-deployment DB state, so the check stays. Documented in the API schema,
`import_vcf --genome-build` and the test-data README; the lab client Phase 4 unblocks should follow it.

## Still open from Phase 3

The grid half of #1706 — top-level specimen and extraction grids under the patients menu, an extraction
link on the sample grid, and private#2837's IDs in the variant page's sample table — is Phase 6, where
the column pass already is. Phase 3 gave both models the URL those links need.

A deployment whose patient CSV ever populated the old Germline/Somatic column gets a manual task out of
`patients/0016` asking a human to decide which of those specimens really were Reference. Nothing
registers on deployments that never used the column.

## Still open from Phase 4

**The client work is raised but not done** —
[`SACGF/variantgrid_api#20`](https://github.com/SACGF/variantgrid_api/issues/20) asks for the
dataclasses and methods behind the patient / specimen / extraction API, the link call, an upload
metadata kwarg (which also reaches Phase 2's `genome_build` and `source`, unreachable from the client
today) and the measures. Deliberately deferred to the end and batched with the other client issues still
to be raised, so it ships as one release rather than several.

The §"Done when" run wants doing once against a real deployment: patient, specimen and both extractions
(`2600000001C` DNA, `2600000001B` RNA) created over the API, the DNA arm's three VCFs reaching
`Sample.extraction` through one seqauto link call, a hand-uploaded file reaching it through upload
metadata alone, a VCF uploaded before its extraction exists parking and then attaching itself, and
TMB/MSI/GIS showing on the specimen page. PR #1716's evidence is the unit suites, which stop short of a
full import with annotation — so this joins Phase 2's #1711 checklist as one pass over a real GRCh37
deployment. Every call is plain JSON over an API token, so the run can be driven by hand rather than
waiting on the client release.

`URLS_NAME_REGISTER` does not gate router URLs, so Shariant still serves the patient API. The names are
in `shariantcommon.py` as intended, but only `api_specimen_measure_bulk_create` is registered through
`perm_path` and therefore actually enforced. Worth a follow-up if Shariant must genuinely not serve
these.

The sample-name fallback cannot name a parent. `Extraction.reference_id` is unique per specimen rather
than globally, so a derived reference matching rows under two specimens parks as Needs attention. A TSO
500 reference embeds its specimen (`2600000001C` starts with `2600000001`), so it does not arise here — a
deployment whose naming does not carry the specimen would need the regex to yield a parent too.

## Phase 5 — gene fusions as variants and the `AllFusions.csv` parser (#1506, phase 1 only)

Independent of Phases 1-4, all of which have landed, so it starts whenever someone is free.

Designed in detail in [`fusion_variants_issue_1506_plan.md`](fusion_variants_issue_1506_plan.md), which
is the spec. In short: gene-level events get a real `Variant`, anchored on a gene-level contig shared by
every build, with the gene pair encoded in a symbolic alt (`<FUSION:HGNC:nnn>`) and a `GeneFusion`
companion model. A `Variant` rather than the bare `Allele` the issue originally proposed, because
`VariantGeneOverlap` and `CohortGenotype` are both `Variant`-keyed — without one, fusions cannot reach
gene lists, compound-het detection or any analysis node. VEP never sees them; a third annotation pipeline
type resolves HGNC → gene locally and writes the overlap rows.

**Scope it to ingestion, storage and identity. Leave fusion classification equivalence alone.** The
user-group feedback on #1506 — "very difficult to implement, needs more thinking and likely more
resourcing, likely a research level project" — is about **fusion equivalence and discordance**, not about
ingesting and displaying them. Somatic classifications already land on
`MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED`, so phase 1 gets that for free and the research-level part
stays deferred.

Two things about the test file that the issue's original design assumed otherwise. Every row carries both
breakpoints (`chr8:128806980` form), so fusions are *not* coordinate-free — what is deferred is the
breakend representation, not the data, which is captured from day one. And the breakpoints are
per-observation rather than per-fusion: `ENTPD3-RPL14` appears three times from one caller with three
different 5′ breakpoints, so identity stays the gene pair and the coordinates live in
`CohortGenotype.info` alongside the rest of the per-row data.

Consequences for the rest of this plan: the fusion file cannot go through the bcftools stages, so the
loader inserts variants directly; and since a `Sample` belongs to exactly one `VCF`, the file creates its
own VCF and Sample rather than joining the RNA arm's `SpliceVariants` sample — Phase 4 ties that sample to
the extraction, Phase 7's grouping node shows both arms together.

**Done when** all 33 rows of the test file import, `EGFR-SEPTIN14` resolves despite the file saying
`SEPT14`, the same gene pair from both callers resolves to one `GeneFusion`, and a gene list containing
`ROS1` finds `CD74-ROS1`.

## Phase 6 — showing non-variants on the grids (#1558, and #1706's grid half)

Less outstanding than the issue title suggests. Per its own comment table, SV already works in the grid
and CNV works as SV — what is missing is surfacing `CohortGenotype.info["CN"]` (and TSO 500's `FORMAT/SM`
linear copy ratio, which importer v21+ already keeps in the format JSON blob but which is not queryable).
That part is independent and could be pulled forward any time.

`cnv.vcf`'s `SEGID` gene symbol joins it, down from Phase 2. The value is stored and the JSON is
queryable, but it is whatever the caller wrote — `MYCL1` where the rest of the pipeline says `MYCL` —
and JSON cannot be joined to `GeneSymbol`/`GeneSymbolAlias`. So it is the same read-time-versus-column
decision as `SM` and `CN`, and wants deciding once for all three.
`GeneSymbolMatcher.get_gene_symbol_id_and_alias_id` (`genes/gene_matching.py:45`) is the resolver
either way.

Fusions are the genuinely new case, and only exist to display after Phase 5. **One grid, not several** —
a fusion is a row like any other, so compound-het detection, gene lists and every downstream node keep
working over the same result set rather than needing a per-kind union. That is also what Phase 5's
choice of a real `Variant` over a bare `Allele` was for. Phase 5 leaves them grid-ready rather than
grid-complete: they carry a gene symbol from their own annotation run, so what remains here is the kind
filter and the columns that only make sense for a fusion.

#1706's grid half joins here because it is the same kind of change and wants the same pass over the
column definitions: top-level specimen and extraction grids under the patients menu, an extraction link
on the sample grid, and private#2837's specimen and patient IDs in the variant page's sample table —
"which of these rows are the same patient" being the question that one is actually asking. Phase 3 gave
both models the URL all of it links to.

## Phase 7 — analysis grouping node (variantgrid_private#223)

Consumes Phase 1 directly. The grouping level is settled: **Extraction, then Specimen** — Patient is too
coarse (same patient sequenced at multiple timepoints, wanted independently) and SequencingSample too
fine (one library/index, so it cannot unite the DNA and RNA arms).

Build the **extraction level first** — same nucleic acid, different callers. That is the VarDict/Mutect
case from sapath#301 and the DRAGEN snv+cnv split, and it has no genome-build complication. Layer
specimen on top after.

Designed in [`tso500_phase_7_plan.md`](tso500_phase_7_plan.md), which is the spec. In short: the node
resolves its samples **at query time** and ORs a per-sample subquery on `pk`, the way `MergeNode`
already does (`analysis/models/nodes/filters/merge_node.py:91-94`), rather than pre-building a
cross-VCF `Cohort`. `Cohort.vcf` being nullable made the cohort look like the obvious vehicle, but
`cohort_genotype_task` hardcodes `format` and `info` to empty and never writes `filters`
(`snpdb/tasks/cohort_genotype_tasks.py:167-173`), so the copy would drop `INFO/CN`, `SEGID`,
`FORMAT/SM` and the FILTER values Phase 2 exists to preserve and Phase 6 exists to display — and its
insert scans all of `snpdb_variant` per extraction. Query-time resolution also has no version to go
stale as Phase 4's out-of-order arrivals set `Sample.extraction`.

It is one node rather than four: `SampleNode` gains a `source_level` and extraction/specimen/patient
FKs, because the only thing that differs between levels is object → sample set, and a separate node
class would force a duplicate of every analysis template. Preserve caller/VCF identity rather than
unioning: sapath#301 requires different cutoffs per caller, and a node that merges everything solves
neither that nor the FFPE-vs-germline case. `analysis_templates_tag`
(`analysis/templatetags/related_analyses_tags.py:134`) takes exactly one of
sample/cohort/trio/quad/pedigree today, with extraction and specimen the natural extension.

Dropping the cohort drops its `genome_build` FK with it, so the single-build constraint that ruled
Patient out is gone — the node restricts to the analysis build, and Patient becomes one more level a
deployment can leave switched off.

## Phase 8 — reporting (sapath#246, #1419, then #444)

Last, because it consumes everything above. #431 notes #246 is old and may be worth rolling into #444.

Designed in [`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md),
which is the spec. In short: somatic reporting sees the same variants over and over, so the phase is
about reusing prior curation. An audit of `EvidenceKey.copy_consensus` first — 142 keys carry it, set
for germline ACMG work and never reviewed against somatic, and several of them copy one patient's
tumour measurements or report narrative onto another's. Then a wizard that triages a run's
`SomaticReportable` tags against what has already been curated for each allele, and walks them one
variant at a time. Then gene-level reuse (#1419), where the human picks the source record because AMP
tiering is gene *and* tumour type and the phenotype data cannot make that call.

- **sapath#246** — keep tagging `SomaticReportable` as now, but launch the classification as AMP/somatic
  from the grid's tags button rather than exporting CSV. Small and independent of the rest; could be
  pulled forward if it is blocking the diagnostic team.
- **#1419** — gene-level copy consensus, deliberately as a copy rather than the first-class gene/disease
  object the issue also proposes. The object needs the disease axis settled first; the copy reuses
  machinery that already works and owns exactly the key set the object would.
- **#444** — multi-variant classification and reporting. Grouped by gene, launched from the sample or
  specimen page, pulling in Phase 4's TMB/MSI measures. `ClassificationReportTemplate` exists; the vue
  template needs to take an array. Starts from the wizard's triaged list.

---

## Parallelism

With Phases 1-4 landed, the model spine is done and what remains splits cleanly: one takes Phase 5
(fusions, `classification` and `upload`), the other Phase 7's grouping node (`analysis`, reading Phase
4's `Sample.extraction`), and they meet at Phase 6, which wants Phase 5's rows in hand. Phase 8 consumes
everything and comes last either way.

With one person, the order above is the order.

## Decisions worth settling before their phase starts

| Decision | Phase | Why it is cheaper now |
|---|---|---|
| Multi-gene partner with one HGNC and one clone identifier — how non-HGNC genes get an identity | 5 | Determines whether `GeneFusion` identity can be non-null on both sides. Parking the row is ruled out: these partners have to be storable |

Phase 6's single-vs-multiple grid question is settled — one grid, so fusions reach comp-het and every
downstream node. Phase 7's four are settled in its own doc.

## Deferred, deliberately

- **Per-gene LOH as a specimen measure** — #1559 lists it, but it is per-gene rather than per-specimen and
  the vendor emits it in a VCF, so it belongs with `abcn_annotated.vcf` below rather than in Phase 4.
- **Fusion equivalence and discordance** — the research-level part of #1506. Somatic classifications get
  `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED` today, so nothing regresses by waiting.
- **Gene/disease curation as a first-class object** — #1419's schema option, and the only thing that can
  express "reviewed on this date, due for review". Needs the disease axis settled; Phase 8 does the
  reuse as a copy in the meantime, over exactly the key set such an object would own.
- **Breakend representation and BND VCF export** — #1506 phases 2 and 3. The breakpoint values themselves
  arrive in `AllFusions.csv` and Phase 5 stores them, so this is a read-side change when a consumer appears.
- **`abcn_annotated.vcf` / gene-level LOH** — cannot start until TAU supplies a file; the format is
  undocumented and not usefully mockable.
- **`_DragenExonCNV.vcf` exact field spelling** — provisional until a run with a real large rearrangement
  arrives. Note the CombinedVariantOutput reports these as `<LOSS>` where the VCF header declares `<DEL>`.
- **`CombinedVariantOutput.tsv`** — settled as not worth ingesting. It is exactly a filtered subset of the
  individual files (`PASS` small variants, `PASS` non-reference CNV) and discards 148 of 149 fusion calls
  and every splice call.
- **Fold-change → DEL/DUP conversion** (sapath#304's open question) — moot for v2.6.2, which emits
  `<DUP>`/`<DEL>` directly with `SM` as the linear copy ratio. The `cnv_tsv_to_vcf.py` command in the
  sapath repo predates that.
