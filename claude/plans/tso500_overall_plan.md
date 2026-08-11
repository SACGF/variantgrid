# TSO 500 — order of work

Sequencing plan for [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431)
and the issues under it. Covers what to build in what order and why, not how to build each piece —
each issue carries its own design.

The file-level ingestion analysis is in [`tso500_ingestion_plan.md`](tso500_ingestion_plan.md); the
test data and its gotchas are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).
Three phases have their own design docs: [`tso500_phase2_plan.md`](tso500_phase2_plan.md) for the loader
refinements, [`fusion_variants_issue_1506_plan.md`](fusion_variants_issue_1506_plan.md) for gene
fusions, and
[`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md) for Phase 8's
reporting.

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
  `issue_1711_tso500_phase2_loader_refinements`, issue #1711). All four items, designed in
  [`tso500_phase2_plan.md`](tso500_phase2_plan.md):

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
  P3  private#2447 tissue status         │                          │
      #1706 pages, search, preview       │                          │
      │              │                   │                          │
  P4  seqauto link ──┤                   │                          │
      #1707 patient API ──► matching ──► #1559 measures ◄───────────┘
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
- `MetricsOutput.tsv` / the run-level TMB summary — feeds Phase 4.
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

## Phase 3 — finish the specimen model, and give it somewhere to show (private#2447, #1706)

Both are Phase 1 leftovers rather than new work, and both are markedly cheaper before Phase 4 puts
`Specimen` behind a public API than after.

### private#2447 — `Specimen.mutation_type` becomes a tissue status

Phase 1 moved `nucleic_acid_source` off `Specimen` onto `Extraction` and left `mutation_type` alone. It
is still a two-value Germline/Somatic field with `default=Mutation.GERMLINE` (`patients/models.py:321`),
and its one substantive consumer is classification autopopulate
(`classification/autopopulate_evidence_keys/evidence_from_sample_and_patient.py:91`), which sets
`SpecialEKeys.ALLELE_ORIGIN` and so feeds `AlleleOriginBucket` and the cross-lab overlaps.

That is a live problem for this work specifically. Every TSO 500 specimen is a tumour block, and unless
somebody sets the field its classifications are stamped "Germline" by default and filed alongside other
labs' real germline records — corrupting exactly the somatic classifications Phase 8 exists to produce.

Two reasons it belongs here rather than "some time":

- Same model, same app, same migration surface as Phase 1, while the context is hot.
- Phase 4 exposes `Specimen` over an API a separate client codes against. Changing the field before that
  is a migration; changing it after is a client contract change.

The design is settled on the issue: a per-specimen `tissue_status` (Reference (unaffected) / Affected
(lesional) / Unknown, defaulting to Unknown), leaving the call-set question to the existing
`VCF.variants_type` and the per-variant question to the existing `allele_origin`, and having autopopulate
assert an origin only where both levels agree. Matched normal stays derivable — "reference specimen, same
patient" — rather than becoming an FK. The migration maps `S` → affected; `G` is ambiguous because it was
the default, so check the deployment's distinct `PatientRecord.specimen_mutation_type` values
(`patients/models.py:662`) before deciding it means anything. The patient CSV's `SPECIMEN_MUTATION_TYPE`
column follows.

### #1706 — specimen and extraction pages, search, preview

`Specimen` and `Extraction` are plain `ExternallyManagedModel` today (`patients/models.py:311`, `:363`):
no `get_absolute_url`, no `PreviewModelMixin`, no search handler, no page of their own. The patient page
tabs (`patients/urls.py:25-26`) are the only place either surfaces.

Phase 4's done-when is "the measures show on the specimen page", so a specimen page is a prerequisite
rather than a nicety. The minimum here: detail pages for both, `PreviewModelMixin` so they hover-card the
way `Patient` does, a search handler, and links each way along `Patient → Specimen → Extraction → Sample`.

The rest of #1706 — top-level specimen/extraction grids under the patients menu, links from the sample
grid, and private#2837's ask for specimen ID on the variant page — is grid work, so it rides with Phase 6
where the grid changes already are. Splitting on that line keeps this phase to what Phase 4 needs.

Beyond #2447's migration this half is views and templates, so it parallelises with anything.

**Done when** a specimen has its own page reachable from search and from the patient, and `tissue_status`
has replaced `mutation_type` everywhere including the CSV.

## Phase 4 — how identifiers cross the lab boundary (#1707, seqauto link, #1559)

One subject in four steps: creating the lab's patients, specimens and extractions in VG, giving an upload
a way to name the extraction it belongs to, absorbing the two arriving in either order, and then hanging
measures off the result.

#431's decision governs the whole phase: **upload the files separately and join later**, not concatenate
before import. So this is assignment, not a new bulk-import path. Concatenating destroys per-caller
cutoffs (which sapath#301 requires), the FORMAT fields across the DRAGEN outputs are disjoint so any merge
is lossy, DNA and RNA are different libraries so a shared sample column makes VAF/depth meaningless, and
fusions are not VCF at all.

There is no run loader to build, either. Upload is strictly per-file and has no directory concept, while
seqauto already models runs, sample sheets and backend VCFs — and the client posts those before it posts a
file. Everything below is assignment on top of the per-file upload that exists.

Raise the client-side issue in `SACGF/variantgrid_api` once the phase lands.

### #1707 — patient / specimen / extraction API

First, because everything after it quotes the identifiers this creates.

`patients` has no serializers at all today; every route into a `Specimen` is a human one (the sample form,
the patient CSV, the admin). The API creates `Patient`, `Specimen` and `Extraction` keyed on the tracking
system's external identifiers, so the client can accession before or alongside posting a run — from a
client that, unlike seqauto, does know the patient.

Do this after Phase 3's `tissue_status` change, so the API ships the field it is keeping.

### Two ways for a VCF to name its extraction

By the time a file arrives the client has already made its seqauto calls (sequencing run, sample sheet,
backend VCFs) and the patient calls above. The file carries a `path`, and on a seqauto deployment that
path is what ties it to the sequencing side: `create_backend_vcf_links` (`upload/vcf/vcf_import.py:409`)
resolves it to a `BackendVCF`, then `link_samples_and_vcfs_to_sequencing` (`:436`) walks that to the
sample sheet and prefix-matches VCF sample names against `SequencingSample` rows
(`get_samples_by_sequencing_sample`, `seqauto/models/models_seqauto.py:1008`), creating
`SampleFromSequencingSample`. So `Sample` ↔ `SequencingSample` already links itself, and the only thing
missing on that path is which extraction the `SequencingSample` is.

**Route 1 — link the `SequencingSample` to the `Extraction`, once per sequencing sample.** Phase 1 added
`SequencingSample.extraction` and nothing sets it from the API. A small link endpoint does, taking the
lookup shape `SequencingSampleLookupSerializer` already defines
(`seqauto/serializers/sequencing_serializers.py:206` — sequencing run, sheet hash, sample name) plus the
extraction identifier. `link_samples_and_vcfs_to_sequencing` then carries it down to `Sample.extraction`
beside the `SampleFromSequencingSample` it already creates, so all three of one arm's VCFs inherit the
extraction from that single call.

Its own endpoint rather than fields on `SampleSheetSerializer`: most sample sheets have no specimen or
extraction to name, and most samples that have an extraction never had a sample sheet, so the two belong
on separate calls. A re-sent sheet builds fresh `SequencingSample` rows (`_create_sequencing_samples`,
`seqauto/serializers/sequencing_serializers.py:250`), so carry an existing extraction across from the
previous current sheet's row with the same `sample_id` rather than making the client re-link.

**Route 2 — an `extraction` key in `FileUpload.metadata`.** Phase 2 built the blob and the per-file-type
key vocabulary, so this is one more key, and it sets `Sample.extraction` with no seqauto records involved
at all — which is what a deployment not running seqauto has, and what a hand-uploaded file has anywhere.
A single-sample VCF takes a bare reference; a multi-sample VCF takes a map keyed on VCF sample name.

Route 1 saves repeating the reference on every file wherever a sample sheet exists; route 2 is the general
mechanism and the only one available off the seqauto path. Where a file has both and they disagree, fail
the import — the same rule Phase 2 applies to a declared `genome_build` the header contradicts.

Validate the key's *shape* at upload time and resolve the reference at import. Phase 2 validates
`genome_build` while the client is still connected precisely so that a typo is a 400; existence-checking
an extraction there would instead reject the ordering race the next section exists to absorb.

### Resolving the identifier — one helper for all three models

The link call, the metadata key and #1707's own payloads all quote the same kind of identifier, and
`Patient`, `Specimen` and `Extraction` are all `ExternallyManagedModel`, so this wants writing once.

An externally managed object is found by its `ExternalPK` (`patients/models.py:65`), which is unique on
`(code, external_type, external_manager)`; one that is not falls back to the local reference Phase 1 put
beside it — `Extraction.reference_id`, `Specimen.reference_id`, `Patient.patient_code`. The resolver
builds its query from whichever of those the caller supplied, rather than insisting on one canonical form:
an external manager and code, an `external_type`, a local reference, or any combination. Each narrows the
query where present and is left out where not.

Nothing needs pinning to a VG-side vocabulary. `external_pk` is a OneToOne *from* each model, so
`Extraction.objects.filter(external_pk__code=...)` can only return extractions — the model being queried
already scopes the lookup, and `external_type` is just one more optional narrowing. The only per-model
knowledge the helper needs is which field holds the local reference, since `Patient` calls it
`patient_code`.

Given both, they narrow to one object and the payload is simply better identified. Given one, it resolves
on that alone. Where the two resolve to *different* objects, that is ambiguity rather than precedence, so
it parks as Needs attention along with everything else the pass cannot settle.

`reference_id` is unique per parent, not globally (`("specimen", "reference_id")`, `("patient",
"reference_id")`), so a bare local reference can legitimately match more than one row. That falls out the
same way — resolve, and park anything ambiguous for a human rather than guessing or making the client send
the whole parent chain. In practice a TSO 500 reference like `2600000001C` already carries its specimen.

Payload shape: a bare string is the code, and an object `{code, external_manager}` disambiguates where a
deployment needs it. The metadata blob is JSON either way, so a multi-sample VCF's per-sample map nests
whichever of the two shapes it needs.

**Keep it on the API rather than the sample-sheet signal.** The SA Pathology app hooks
`sequencing_run_sample_sheet_created_signal` (`seqauto/models/models_seqauto.py:248`) to assign patients
from HelixIDs in the sheet, and the same hook could reach extractions. Leave it where it is: it only fires
where seqauto runs, it is per-deployment glue rather than anything an API contract can state, and route 2
has to exist regardless.

### Reconciling what arrives out of order — Mocha's pattern

Mocha hit this and the answer is worth copying wholesale (SACGF/mocha#126, commits `ed1d0df`, `982022e`):
independent feeds that legitimately arrive in any order, and callers who should never have to re-send.
What it does:

- **Nullable FK, with `match_status` (Matched / Pending / Needs attention) and a human-readable
  `match_error` on the same row.** A row that cannot resolve yet is parked rather than rejected.
- **Reject only what can never resolve** — no anchor record for the identifiers at all — with a 400.
  Everything else is a 202 carrying the reason, so the caller knows it is pending rather than failed.
- **Key idempotency on the posted identifiers, never on the nullable FK.** Postgres treats nulls as
  distinct, so a `unique_together` spanning the FK lets a parked resend duplicate.
- **A periodic `reconcile_pending_*` task** re-resolves parked rows, fired on schedule and again whenever
  new upstream data lands. Its `apply_match` leaves an already-matched row alone, so a confirmed link
  never flaps back.
- **Promote Pending to Needs attention after a few days** — past that window it is a real mismatch rather
  than the load race, and wants a human instead of rotting silently.
- **Surface it**: `status` and `match_error` read-only on the serializer, and unmatched rows filterable on
  the list page.

Here the parked claim always has a row to live on, so it needs no new table: an unresolved reference sits
beside the `extraction` FK that already exists on both `Sample` and `SequencingSample`, as
`extraction_reference` plus the two match columns. Reconciling is then one task across both models, over
the rows whose FK is null and whose reference is set.

Where no tracking system supplies identifiers at all, parsing the sample name is the fallback —
`2600000001` is the specimen and the trailing `C`/`B` name the DNA and RNA extractions. A fallback, not
the mechanism: it is the only option for those deployments, and guesswork everywhere else.

### #1559 — specimen-level measures

Small, and it proves Phase 1's identifiers work across the LIMS boundary. Direction is already settled:
**these arrive by API from the lab client, not by VG parsing pipeline output.** The measures are
scattered across the run — TMB summary in one place, MSI in a separate JSON, GIS with tumour fraction
and ploidy in another written by a different tool — and none of that structure is worth VG knowing about.
When Illumina moves a file between releases it becomes a client change, not a VG release.

Four things to hold to, from the issue:

- **The client transcribes, it never computes.** Lifting a value out of vendor output is fine; deriving
  one moves a clinical calculation into an unversioned client nobody reviews.
- **Store the raw payload beside the parsed fields**, so "which file and which pipeline version produced
  this number" is answerable at report time. `HelixNGSOrder.raw_data` in Mocha is the existing pattern.
- **Key on the specimen/extraction external identifiers from #1704**, not a vendor pair-ID string.
- **Keep the schema vendor-neutral** — `SpecimenMeasure(specimen, measure_type, value, unit, method,
  source_payload)` rather than a TSO500-shaped table.

Store the **score as well as the call**. HRD gives a genomic instability score and no positive/negative
call — the threshold is lab policy, not vendor output. Same for MSI and TMB high/low. Without the raw
value plus the threshold applied and by whom, a later re-interpretation is unreconstructable. This
mirrors the evidence keys already in `classification/migrations/0164_somatic_hrd_msi_tmb_ekeys.py` —
`somatic:tmb_value` beside `somatic:tmb_status`, `somatic:msi_value` beside `somatic:msi_status`.
`somatic:hrd_status` has no value counterpart yet; a score field is the equivalent.

Per-gene absolute and minor copy number are *not* in scope here — they are per-gene rather than
per-specimen and arrive as a VCF, so they are ingestion. They land whenever Phase 0's request produces
`abcn_annotated.vcf`.

**Done when** the whole test run lands against one specimen with both arms attached — patient, specimen
and extractions created over the API, the DNA arm's VCFs reaching `Sample.extraction` through one seqauto
link call and a hand-uploaded file reaching it through upload metadata alone — and TMB/MSI/GIS posted
against that specimen shows on its page.

## Phase 5 — gene fusions as variants and the `AllFusions.csv` parser (#1506, phase 1 only)

Independent of Phases 1-4, so this can run in parallel with a second person from Phase 0 onward.

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

Fusions are the genuinely new case, and only exist to display after Phase 5. The single-grid vs
multiple-grids question is worth deciding on real fusion rows rather than in the abstract — with fusion
`Variant`s in hand, "sort by gene brings the fusions together" is testable instead of hypothetical. Phase
5 leaves them grid-ready rather than grid-complete: they carry a gene symbol from their own annotation
run, so what remains here is the kind filter and the columns that only make sense for a fusion.

#1706's grid half joins here because it is the same kind of change and wants the same pass over the
column definitions: top-level specimen and extraction grids under the patients menu, an extraction link
on the sample grid, and private#2837's specimen and patient IDs in the variant page's sample table —
"which of these rows are the same patient" being the question that one is actually asking. All of it
depends on Phase 3 having given both models a URL to link to.

## Phase 7 — analysis grouping node (variantgrid_private#223)

Consumes Phase 1 directly. The grouping level is settled: **Extraction, then Specimen** — Patient is too
coarse (same patient sequenced at multiple timepoints, wanted independently) and SequencingSample too
fine (one library/index, so it cannot unite the DNA and RNA arms).

Build the **extraction level first** — same nucleic acid, different callers. That is the VarDict/Mutect
case from sapath#301 and the DRAGEN snv+cnv split, and it has no genome-build complication. Layer
specimen on top after.

Likely less new machinery than it looks: `Cohort.vcf` is already nullable and spans VCFs
(`snpdb/models/models_cohort.py:53`), `CohortNode` consumes it today, and the per-caller-cutoff case is
already expressible as N× `SampleNode` with their own `min_ad`/`min_dp` into a `MergeNode`. So this is
largely an entry-point problem — auto-selecting the right samples and wiring the template — and
`analysis_templates_tag` (`analysis/templatetags/related_analyses_tags.py:134`) takes exactly one of
sample/cohort/trio/quad/pedigree today, with extraction and specimen the natural extension. Preserve
caller/VCF identity rather than unioning: sapath#301 requires different cutoffs per caller, and a node
that merges everything solves neither that nor the FFPE-vs-germline case.

`Cohort` has a `genome_build` FK, so a cross-VCF cohort is single-build by construction — fine at
extraction and specimen level, and only a problem if this ever goes up to Patient.

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

With two people: one takes Phase 3 → 4 (the model spine, now that Phase 1 has landed), the other takes
Phase 5 (fusions), and they meet at Phase 6. Phase 5 touches `classification` and `upload`; Phases 3-4
touch `patients`, `snpdb` and `seqauto`, so the conflict surface is small.

Phase 3's #1706 half is views and templates only, so it splits off cleanly on top of that. One ordering
constraint remains across the lot: #2447's field change lands before #1707 serializes it. Phase 4's
route 2 is unblocked either way — the `FileUpload.metadata` blob it hangs an `extraction` key off
already exists, along with the per-file-type key declaration and upload-time validation it needs.

With one person, the order above is the order.

## Decisions worth settling before their phase starts

| Decision | Phase | Why it is cheaper now |
|---|---|---|
| Does `G` in `PatientRecord.specimen_mutation_type` mean anything on each deployment? | 3 | Decides whether the migration can map `G` → Unknown and lose nothing |
| Multi-gene partner with one HGNC and one clone identifier — park the row or resolve to the HGNC member? | 5 | Determines whether `GeneFusion` identity can be non-null on both sides |
| Single grid or multiple grids for non-variants? | 6 | Best answered against real fusion rows, so genuinely wait for Phase 5 |

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
