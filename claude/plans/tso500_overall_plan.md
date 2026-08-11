# TSO 500 — order of work

Sequencing plan for [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431)
and the issues under it. Covers what to build in what order and why, not how to build each piece —
each issue carries its own design.

The file-level ingestion analysis is in [`tso500_ingestion_plan.md`](tso500_ingestion_plan.md); the
test data and its gotchas are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).

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

## Where things stand

The five test files from `ed5e15a33` — one de-identified run, one specimen, DNA and RNA arms:

| File | State |
|---|---|
| `hard-filtered.vcf` | loads, 93 records |
| `cnv.vcf` | loads, 25 records — copy-neutral rows and `SM` are Phase 2 refinements |
| `SpliceVariants.vcf` | loads, 17 records — `AD`/`DP` mapping and the dropped filter are Phase 2 |
| `_DragenExonCNV.vcf` | reaches the importer after Phase 0; still needs an explicit genome build (Phase 2) |
| `AllFusions.csv` | no longer breaks file-type detection; needs a parser and somewhere to put the rows (Phase 5) |

## Dependency map

```
  Phase 0  pipeline fixes ✔ ──► P2  loader refinements   (no model deps — anyone, any time)
           chase TAU files ─────────────────────────────────────────┐
                                                                    │
  Phase 1  #1704 Specimen→Extraction ✔                              │
             │                                                      │
      ┌──────┴───────────────────────────┐                          │
      │                                  │                          │
  P3  private#2447 tissue status         │                          │
      #1706 pages, search, preview       │                          │
      │              │                   │                          │
  P4  seqauto text ──┤                   │                          │
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

## Phase 2 — per-file loader refinements

Make each file import with the right values in it. Nothing here depends on Phase 1, or on anything else
in this plan — four independent items, each testable against one test file. It could have gone in
Phase 0 and it can be picked up by anyone at any point, which is why it is separated from the linkage
work that used to share this phase (now Phase 4).

Designed in detail in [`tso500_phase2_plan.md`](tso500_phase2_plan.md), which is the spec.

- **Upload metadata — genome build and source.** `_DragenExonCNV.vcf` has no `##contig` lines and only
  `##reference=file:resource_bundle/hg19_decoy/`, so `vcf_detect_genome_build_from_header` finds nothing
  and `upload/tasks/vcf/genotype_vcf_tasks.py:63` sets `REQUIRES_USER_INPUT` and stops; `cnv.vcf` and
  `hard-filtered.vcf` declare no `##source`, so nothing keyed on the caller can reach them. Whatever
  submits the file knows both, so let it say so explicitly instead of relying on header detection.
  A `FileUpload.metadata` JSON blob, validated at upload time.
- **`VCFSourceSettings` mapping for SpliceGirl.** The `##source` line (`SpliceGirl 1.0.0.614`) is already
  the hook. In this file `AD` is `Number=1` splice-supporting reads and `DP` is *reference* reads, so
  binding them by name silently gives empty allele depth, missing VAF, and a read-depth column showing
  reference counts. Derive VAF from `ALTDEDUP`/`REFDEDUP`.
- **Preserve `LowUniqueAlignments`.** The splice header declares `LowQ` only, so `vcf_clean_and_filter`
  strips the undeclared code from 15 of 17 records — and that is the filter separating the two `PASS`
  oncology calls from the background. The genotype processor already copes via `_add_undeclared_filter_code`,
  so the fix is letting undeclared values through the header stage.
- **Skip copy-neutral `cnv.vcf` rows.** `ALT=.` is 406 of 500 in a real run; they carry no call and
  currently import as reference variants with the segment span discarded.
Resolving `cnv.vcf`'s `SEGID` gene symbols through `GeneSymbol` was a fifth item here, and has moved to
Phase 6: `CohortGenotype.info` is a `JSONField`, so nothing is lost at import, and whether to resolve
aliases at read time or into a queryable column is the same question `SM` and `CN` pose there.

None of these have their own issue — they came out of running the test data through the real pipe stages,
and sapath#431 is the only thing above them. Worth raising individually if they are going to be picked up
by different people.

**Done when** all five test files import with the right values: splice VAF derived from
`ALTDEDUP`/`REFDEDUP`, `LowUniqueAlignments` surviving to the grid, no copy-neutral rows inserted, and
`_DragenExonCNV.vcf` landing against a declared build.

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

## Phase 4 — how identifiers cross the lab boundary (seqauto linkage, #1707, #1559)

One subject in four steps: getting the lab's specimen and extraction identifiers into VG, turning them
into objects, reconciling the two, and then hanging measures off the result. The sample-sheet half used to
sit in Phase 2, but it is the same boundary and the same client as #1707, and it reads as one design only
when the two are together.

#431's decision governs the whole phase: **upload the files separately and join later**, not concatenate
before import. So this is assignment, not a new bulk-import path. Concatenating destroys per-caller
cutoffs (which sapath#301 requires), the FORMAT fields across the DRAGEN outputs are disjoint so any merge
is lossy, DNA and RNA are different libraries so a shared sample column makes VAF/depth meaningless, and
fusions are not VCF at all.

Raise the client-side issue in `SACGF/variantgrid_api` once the phase lands.

### The sample sheet — carry `Specimen`/`Extraction` as text beside the FKs

Phase 1 added `SequencingSample.extraction`, but `SequencingSampleSerializer`
(`seqauto/serializers/sequencing_serializers.py:232`) does not list the field, so the FK is unreachable
from the lab client that posts sample sheets — the only ways to set an extraction today are the sample
form, the patient CSV and the admin.

**Settled**: `SequencingSample` grows plain-text specimen and extraction reference columns alongside the
model FKs. The client sends the text, and the API stores whatever it was sent — the payload is always
accepted. Where a `Specimen`/`Extraction` already exists the serializer links it; where one does not, the
text is kept and the FK stays null. The lab client owns the identifiers, VG owns the objects, and a run
posts successfully whether or not the tracking system got there first. This is the shape `PatientRecord`
already uses (`patients/models.py:635`) — raw columns from the file beside the matched FK.

Because nothing is created, this step needs no `Patient`, which is what would otherwise have blocked it:
`Specimen.patient` is non-null and nothing in the seqauto payload names a patient. That gap is the next
step's job.

### #1707 — patient / specimen / extraction API

`patients` has no serializers at all today; every route into a `Specimen` is a human one (the sample form,
the patient CSV, the admin). The API creates `Patient`, `Specimen` and `Extraction` keyed on the tracking
system's external identifiers, so the client can accession before or alongside posting a run — from a
client that, unlike seqauto, does know the patient.

Do this after Phase 3's `tissue_status` change, so the API ships the field it is keeping.

### Reconciling the two — the matching pass

With both halves in, the pass has real objects to match against rather than hypothetical ones: stored
sample-sheet text against specimens and extractions that arrived afterwards, on the R5 external
identifiers where the tracking system supplies them.

Where it does not, fall back to parsing the sample name — `2600000001` is the specimen and the trailing
`C`/`B` name the DNA and RNA extractions. A fallback, not the mechanism: it is the only option for
deployments with no tracking system, and it is guesswork everywhere else.

**Decide before starting:** whether anything beyond this is needed to get a run's files loaded. Upload is
strictly per-file today and has no directory concept; seqauto already models runs, sample sheets and file
discovery. If the client posts the linkage and the files go through normal upload, there is no run loader
to site — which is the answer this phase is built around, so it is worth confirming rather than assuming.

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

**Done when** the whole test run lands against one specimen with both arms attached — sample sheet posted,
patient/specimen/extraction created over the API, FKs resolved — and TMB/MSI/GIS posted against that
specimen shows on its page.

## Phase 5 — `GeneFusion` and the `AllFusions.csv` parser (#1506, phase 1 only)

Independent of Phases 1-4, so this can run in parallel with a second person from Phase 0 onward.

**Scope it to ingestion, storage and identity. Leave fusion classification equivalence alone.** The
design on #1506 already does this deliberately: a `GeneFusion` model keyed on a canonical gene pair with
a one-to-one *bare* `Allele` (no `ClinGenAllele`, no `VariantAllele` — purely a grouping key), plus an
`ImportedAlleleInfo.gene_fusion` FK that short-circuits HGVS resolution. Because the bare Allele slots
into the existing infrastructure, `ClinicalContext`, `DiscordanceReport`, the VCF pipeline, annotation
and ClinGen all need no changes; range queries just exclude `genefusion__isnull=False`.

This matters for sequencing, because the user-group feedback on #1506 — "very difficult to implement,
needs more thinking and likely more resourcing, likely a research level project" — is about **fusion
equivalence and discordance**, not about ingesting and displaying them. Somatic classifications already
land on `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED`, so phase 1 gets that for free and the research-level
part stays deferred. Breakend coordinates (#1506 phase 2) and BND VCF export (phase 3) are also deferred.

The parser then needs a fusion import task factory claiming the file — recognisable from the
`# Source = FusionProcessor` first line and the `Caller,Gene A,Gene B,…` header row — with a processing
ability above `GeneListImportTaskFactory`'s 1, which would otherwise win and import the fusions as a
gene list. Real-file properties it has to survive, all present in the test data:

- **Multi-gene partners are routine**, not edge cases — `ROS1;GOPC`, `RP11-458D21.5;NOTCH2NL`. The
  per-gene annotation columns are correspondingly semicolon-delimited. **Open question:** for a partner
  that is one HGNC gene plus one clone-based identifier, reject the row or resolve to the HGNC member?
- **`SEPT14`** (EGFR-SEPT14) — HGNC renamed it `SEPTIN14`, and spreadsheets date-mangle it. A direct test
  of keying on HGNC ID rather than symbol.
- **Mixed delimiters** — `PPARG/AC016683.6` uses a slash alongside the `;` and `-` forms.
- **Direction is reported, not inferred** — `Fusion Directionality Known` plus `Gene A Sense`/`Gene B Sense`
  populate `is_ordered` from data.
- **Breakpoint context** — `Gene A Location`/`Gene B Location` carry `IntactExon`/`BrokenExon`/`Intronic`/
  `Intergenic`. Keep in `source_data`; a later breakend representation wants exactly this and it costs
  nothing now.
- **Two callers in one file** (`DRAGEN`, `SpliceGirl`) whose scores and filters the header warns are not
  comparable — carry caller per row.
- **Ingest unfiltered.** 1 of 149 rows passed the caller's own filter. Same posture as the CNV and
  small-variant VCFs: take everything, apply our own thresholds.

**Done when** all 33 rows of the test file import, `EGFR-SEPTIN14` resolves despite the file saying
`SEPT14`, and the same fusion from both callers resolves to one `GeneFusion`.

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
multiple-grids question is worth deciding on real fusion rows rather than in the abstract — with a
`GeneFusion`-keyed `Allele` in hand, "sort by gene brings the fusions together" is testable instead of
hypothetical.

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

## Phase 8 — reporting (sapath#246, then #444)

Last, because it consumes everything above. #431 notes #246 is old and may be worth rolling into #444.

- **sapath#246** — keep tagging `SomaticReportable` as now, but launch the classification as AMP/somatic
  from the grid's tags button rather than exporting CSV. Small and independent of the rest; could be
  pulled forward if it is blocking the diagnostic team.
- **#444** — multi-variant classification and reporting. Grouped by gene, launched from the sample or
  specimen page, pulling in Phase 4's TMB/MSI measures. `ClassificationReportTemplate` exists; the vue
  template needs to take an array.

---

## Parallelism

With two people: one takes Phase 3 → 4 (the model spine, now that Phase 1 has landed), the other takes
Phase 5 (fusions), and they meet at Phase 6. Phase 5 touches `classification` and `upload`; Phases 3-4
touch `patients`, `snpdb` and `seqauto`, so the conflict surface is small.

Two things split off cleanly on top of that. Phase 2 is four independent file-correctness fixes, so it
can go to whoever, whenever — including before Phase 3. Phase 3's #1706 half is
views and templates only. The one ordering constraint across the lot is that #2447's field change lands
before #1707 serializes it.

With one person, the order above is the order.

## Decisions worth settling before their phase starts

| Decision | Phase | Why it is cheaper now |
|---|---|---|
| Is a run loader needed at all, or does client-posted linkage plus normal upload cover it? | 4 | The whole phase is built around the second answer; confirm rather than assume |
| Does `G` in `PatientRecord.specimen_mutation_type` mean anything on each deployment? | 3 | Decides whether the migration can map `G` → Unknown and lose nothing |
| Multi-gene partner with one HGNC and one clone identifier — reject or resolve? | 5 | Determines whether `GeneFusion` identity can be non-null on both sides |
| Single grid or multiple grids for non-variants? | 6 | Best answered against real fusion rows, so genuinely wait for Phase 5 |

## Deferred, deliberately

- **Per-gene LOH as a specimen measure** — #1559 lists it, but it is per-gene rather than per-specimen and
  the vendor emits it in a VCF, so it belongs with `abcn_annotated.vcf` below rather than in Phase 4.
- **Fusion equivalence and discordance** — the research-level part of #1506. Somatic classifications get
  `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED` today, so nothing regresses by waiting.
- **Breakend coordinates and BND VCF export** — #1506 phases 2 and 3.
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
