# TSO 500 — order of work

Sequencing plan for [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431)
and the issues under it. Covers what to build in what order and why, not how to build each piece —
each issue carries its own design.

The file-level ingestion analysis is in [`tso500_ingestion_plan.md`](tso500_ingestion_plan.md); the
test data and its gotchas are in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).

---

## Where things stand

Test data landed in `ed5e15a33` — one de-identified run, one specimen, DNA and RNA arms, five files.
Running each through the real pipe stages gave:

| File | Today | Blocked on |
|---|---|---|
| `hard-filtered.vcf` | loads, 93 records | — |
| `cnv.vcf` | loads, 25 records | copy-neutral rows and `SM` are refinements |
| `SpliceVariants.vcf` | loads, 17 records | `AD`/`DP` mean something else; filter dropped |
| `_DragenExonCNV.vcf` | dies at first pipe stage | PyVCF FILTER parse; no `##contig` lines |
| `AllFusions.csv` | errors during file-type detection | no parser, and nowhere to put the rows |

So three of five files already load, and the two that don't are blocked by different things — one by a
small pipeline bug, one by a missing model. That split is what sets the order below.

## The single most useful thing to fix first

Everything downstream keys on the same thing: **the identifier pair that says these two arms are one
specimen**. In the test run that is accession `2600000001` with container suffixes `C` (DNA) and `B`
(RNA). Four separate pieces of work need that join to exist:

- the DNA and RNA arms have to be recognisable as one specimen at all (#431's "C- ID / Pair ID")
- specimen-level measures arrive by API keyed on it, not on a vendor pair-ID string
  ([#1559](https://github.com/SACGF/variantgrid/issues/1559))
- the analysis node selects on it
  ([variantgrid_private#223](https://github.com/SACGF/variantgrid_private/issues/223) — settled on
  Extraction/Specimen as the grouping level, since Patient is too coarse and SequencingSample too fine)
- multi-variant reporting groups by it ([#444](https://github.com/SACGF/variantgrid/issues/444))

That join key is [#1704](https://github.com/SACGF/variantgrid/issues/1704). It is the right next step,
and the reason is not just that other issues block on it — it is that `Sample.specimen` is barely
wired in today (five consumers, listed in the issue's R7) and every month of somatic work makes that
surface larger. It is the one change here that gets more expensive by waiting.

## Dependency map

```
  Phase 0  pipeline fixes ─────────────────┐
           chase TAU files ────────────┐   │
                                       │   │
  Phase 1  #1704 Specimen→Extraction   │   │
             │                         │   │
      ┌──────┼───────────────┐         │   │
      │      │               │         │   │
  P2  │  TSO500 linking ◄────┼─────────┼───┘
      │      + loader refinements       │
      │                                 │
  P3  │  #1559 specimen measures (API)  │
      │                                 │
  P4  #1506 GeneFusion ──► AllFusions parser ◄──┘   (independent of #1704 — can run in parallel)
      │
  P5  #1558 grid display
      │
  P6  private#223 analysis grouping node
      │
  P7  sapath#246 ──► #444 multi-variant reporting
```

Only two things are genuinely on a critical path: #1704 before anything that groups samples, and
#1506 before the fusion parser has a destination. Everything else is ordered by value, not by
necessity.

---

## Phase 0 — unblock the remaining test data, and start the long-lead requests

Two small pipeline fixes, both worth doing regardless of TSO 500, both prerequisites for exercising
the test data end to end.

1. **Read declared FILTER IDs off the raw header lines** in `upload/management/commands/vcf_clean_and_filter.py`.
   Illumina's LrCalculator writes `##FILTER` lines with keys beyond `ID`/`Description`, which PyVCF
   rejects outright, killing the whole import at the first pipe stage. This is the last PyVCF use in the pipe.
2. **`get_processing_ability` returns 0 for unreadable input** —
   `upload/import_task_factories/import_task_factories.py:127` calls `pd.read_csv` bare, so one CSV
   shape it cannot parse breaks file-type detection for every other CSV type.

**Done** — `vcf_header_filter_ids()` in `library/genomics/vcf_utils.py`, tests in
`library/tests/test_vcf_header_filter_ids.py`. Note this went to a regex over the header lines rather
than cyvcf2 as originally sketched: cyvcf2 needs a real file descriptor (`StringIO`/`BytesIO` both raise
`UnsupportedOperation: fileno`), and this command is a streaming stdin filter, so spilling the header to
a temp file just to read it back would be worse than the hand-rolled parse the rest of the command
already uses — which exists precisely because strict parsers choke on real-world headers.

In parallel, and starting now because they have lead time, **request the missing files from TAU**:

- `Logs_Intermediates/Gis/<sample>/<sample>.abcn_annotated.vcf` + `<sample>.abcn_genes.tsv` — per-gene
  absolute and minor copy number, the latter being where gene-level LOH comes from. The only genuinely
  CVO-only data, and a VCF, so it is in ingestion scope. Illumina names the path but publishes no
  INFO/FORMAT spec, so it cannot be mocked up.
- `MetricsOutput.tsv` / the run-level TMB summary — feeds Phase 3.
- A run with a real BRCA1/2 large rearrangement. Both `_DragenExonCNV.vcf` records in the test data are
  constructed, and no published TSO 500 output has a populated one either. Until one arrives the field
  spelling in that file stays provisional.

**Done when** all five test files at least reach the importer without erroring, and the requests are out.

## Phase 1 — `Specimen → Extraction → Sample` (#1704)

The foundational model change. Read the issue for R1-R8; the shape is:

```
Patient └── Specimen └── Extraction └── [SequencingSample] └── Sample
```

Points worth holding onto while implementing:

- **`Extraction` goes in `patients`, not `seqauto`** (R1). VCFs arrive by hand, and some deployments
  never run seqauto, so grouping has to work without it.
- **`nucleic_acid_source` moves to `Extraction`** (R3). It sitting on `Specimen` is the direct cause of
  the DNA and RNA arms becoming two unrelated specimens. `mutation_type` stays on `Specimen` — it
  describes the material. `Sample.variants_type` is a third, genuinely different thing and stays put.
- **`Specimen` needs a surrogate PK** (R4). `reference_id` is a `TextField` primary key today, so
  specimen references are globally unique across all patients and any FK drags text keys around.
- **Both levels extend `ExternallyManagedModel`** (`patients/models.py:85`) (R5) — this is what lets VG
  and the LIMS agree on a specimen without string surgery on sample names, and it is what Phase 3 keys on.

Consumer surface, all of it: `patients/import_records.py`, `snpdb/forms.py:366`, `snpdb/grids.py:161`,
`patients/views_autocomplete.py`, and `get_evidence_fields_for_specimen()` in
`classification/autopopulate_evidence_keys/evidence_from_sample_and_patient.py:86`. Patient CSV import
has to keep round-tripping `PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE` (R6) — decide whether the
existing column implies an extraction or the CSV format grows one.

**Settled** (2026-08-10), closing the open questions on the issue:

- **Name it `Extraction`.** The LIMS-side collision noted on #1704 is real but lives in documentation,
  not in the model.
- **No QC fields on `Extraction`** (yield, concentration, 260/280). R5's `external_pk` keeps QC
  reconcilable in the tracking system, and VG never generates those numbers — a second copy would only
  drift. Nullable columns stay cheap to add if a consumer appears.
- **No parent above `Specimen`.** A block (tumour FFPE) *is* a specimen, so several blocks from one
  clinical case are several `Specimen` rows on one `Patient`, already told apart by `tissue` and
  `collection_date`. A case level would only earn its place if one case needed its blocks analysed
  together, which TSO 500 does not.
- **The existing patient CSV column implies a single extraction** — `SPECIMEN_NUCLEIC_ACID_SOURCE`
  creates/updates one `Extraction` under the `Specimen`, so existing CSVs keep working untouched (R6).
  Those files describe one extraction per specimen anyway.

Existing `Specimen` rows need a data migration creating one `Extraction` each, carrying the specimen's
current `nucleic_acid_source`, and repointing their samples onto it.

**Done when** a `Sample` reaches its specimen via `sample.extraction.specimen`, one specimen holds a DNA
and an RNA extraction, and the existing patient CSV import still round-trips.

**Built** on branch `issue_1704_specimen_extraction_sample` — models, the four migrations, R7's consumers,
and an `Extraction` formset on the patient specimens tab so DNA/RNA stays settable now that it has left
the specimen form. Two surfaces stayed specimen-shaped and are worth revisiting as consumers appear: the
patient CSV columns are all `SPECIMEN_*`, so a CSV can name one extraction per specimen but not the DNA
and RNA arms of one block; and `specimen_autocomplete` no longer has a caller, the sample form having
moved to `extraction_autocomplete` (kept, since the sapath and private repos may use it).

## Phase 2 — TSO 500 linking, and the per-file loader refinements

With the model in place, make the run's files land against the right extractions. #431's decision is
**upload the files separately and join later** — not concatenate before import — so this is mostly
assignment, not a new bulk-import path. That is the cheap version and it is the right one: concatenating
destroys per-caller cutoffs (which sapath#301 requires), the FORMAT fields across the DRAGEN outputs are
disjoint so any merge is lossy, DNA and RNA are different libraries so a shared sample column makes
VAF/depth meaningless, and fusions are not VCF at all.

- **Parse accession + container suffix into `Specimen` + `Extraction`.** `2600000001` is the specimen,
  the trailing `C`/`B` name the two extractions. Match on external identifiers where the tracking system
  supplies them rather than trusting the name parse.
- **Supply the genome build explicitly.** `_DragenExonCNV.vcf` has no `##contig` lines and only
  `##reference=file:resource_bundle/hg19_decoy/`, so `vcf_detect_genome_build_from_header` finds nothing
  and `upload/tasks/vcf/genotype_vcf_tasks.py:63` sets `REQUIRES_USER_INPUT` and stops. The loader knows
  the build; pass it alongside the file.
- **`VCFSourceSettings` mapping for SpliceGirl.** The `##source` line (`SpliceGirl 1.0.0.614`) is already
  the hook. In this file `AD` is `Number=1` splice-supporting reads and `DP` is *reference* reads, so
  binding them by name silently gives empty allele depth, missing VAF, and a read-depth column showing
  reference counts. Derive VAF from `ALTDEDUP`/`REFDEDUP`.
- **Preserve `LowUniqueAlignments`.** The splice header declares `LowQ` only, so `vcf_clean_and_filter`
  strips the undeclared code from 15 of 17 records — and that is the filter separating the two `PASS`
  oncology calls from the background. The genotype processor already copes via `_add_undeclared_filter_code`;
  either let undeclared values through the header stage or have the loader add the `##FILTER` line.
- **Skip copy-neutral `cnv.vcf` rows.** `ALT=.` is 406 of 500 in a real run; they carry no call and
  currently import as reference variants with the segment span discarded.
- **Carry `Specimen`/`Extraction` through the seqauto API as text beside the FKs.** Phase 1 added
  `SequencingSample.extraction`, but `SequencingSampleSerializer`
  (`seqauto/serializers/sequencing_serializers.py:232`) does not list the field, so the FK is
  unreachable from the lab client that posts sample sheets — the only ways to set an extraction today
  are the sample form, the patient CSV and the admin.

  **Settled**: `SequencingSample` grows plain-text specimen and extraction reference columns alongside
  the model FKs. The client sends the text, and the API stores whatever it was sent — the payload is
  always accepted. Where a `Specimen`/`Extraction` already exists the serializer links it; where one
  does not, the text is kept and the FK stays null for a later matching pass. The lab client owns the
  identifiers, VG owns the objects, and a run posts successfully whether or not the tracking system got
  there first. This is the shape `PatientRecord` already uses (`patients/models.py:635`) — raw columns
  from the file beside the matched FK — so the matching pass has a precedent to follow.

  Because nothing is created, the API needs no `Patient`, which is what would otherwise have blocked
  this: `Specimen.patient` is non-null and nothing in the seqauto payload names a patient.

  **Deferred to later work**: the matching pass itself — reconciling stored text against specimens and
  extractions that arrive afterwards, on the R5 external identifiers where the tracking system supplies
  them and the accession/container-suffix parse above where it does not.

Note `SEGID=MYCL1` uses an older symbol than the rest of the pipeline (`MYCL`) — resolve gene symbols
through `GeneSymbol` rather than trusting what the caller wrote. The same argument recurs in Phase 4.

**Decide before starting:** where the loader lives. Upload is strictly per-file today and has no
directory concept; seqauto already models runs, sample sheets and file discovery. Given "upload
separately, join later", the lightest answer is probably that the lab client posts the linkage and the
files go through normal upload — which also matches Phase 3's direction. Worth settling explicitly
rather than by default.

**Done when** the whole test run imports with both arms attached to one specimen, splice VAF is right,
and no copy-neutral rows are inserted.

## Phase 3 — specimen-level measures (#1559)

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

**Done when** the client can post TMB/MSI/GIS for a specimen and it shows on the specimen page.

## Phase 4 — `GeneFusion` and the `AllFusions.csv` parser (#1506, phase 1 only)

Independent of Phases 1-3, so this can run in parallel with a second person from Phase 0 onward.

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

## Phase 5 — showing non-variants on the grids (#1558)

Less outstanding than the issue title suggests. Per its own comment table, SV already works in the grid
and CNV works as SV — what is missing is surfacing `CohortGenotype.info["CN"]` (and TSO 500's `FORMAT/SM`
linear copy ratio, which importer v21+ already keeps in the format JSON blob but which is not queryable).
That part is independent and could be pulled forward any time.

Fusions are the genuinely new case, and only exist to display after Phase 4. The single-grid vs
multiple-grids question is worth deciding on real fusion rows rather than in the abstract — with a
`GeneFusion`-keyed `Allele` in hand, "sort by gene brings the fusions together" is testable instead of
hypothetical.

## Phase 6 — analysis grouping node (variantgrid_private#223)

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

## Phase 7 — reporting (sapath#246, then #444)

Last, because it consumes everything above. #431 notes #246 is old and may be worth rolling into #444.

- **sapath#246** — keep tagging `SomaticReportable` as now, but launch the classification as AMP/somatic
  from the grid's tags button rather than exporting CSV. Small and independent of the rest; could be
  pulled forward if it is blocking the diagnostic team.
- **#444** — multi-variant classification and reporting. Grouped by gene, launched from the sample or
  specimen page, pulling in Phase 3's TMB/MSI measures. `ClassificationReportTemplate` exists; the vue
  template needs to take an array.

---

## Parallelism

With two people: one takes Phase 1 → 2 → 3 (the model spine), the other takes Phase 0 → 4 (fusions), and
they meet at Phase 5. Phase 4 touches `classification` and `upload`; Phases 1-3 touch `patients` and
`snpdb`, so the conflict surface is small.

With one person, the order above is the order.

## Decisions worth settling before their phase starts

| Decision | Phase | Why it is cheaper now |
|---|---|---|
| ~~Extraction QC, Specimen parent, patient CSV, naming~~ | 1 | Settled — see Phase 1 |
| Where does the run loader live — upload, seqauto, or lab-client API? | 2 | Settle explicitly; the API answer lines up with Phase 3 |
| ~~May the seqauto API create specimens/extractions, or only match them?~~ | 2 | Settled — text beside the FKs, link what exists, match the rest later. See Phase 2 |
| Multi-gene partner with one HGNC and one clone identifier — reject or resolve? | 4 | Determines whether `GeneFusion` identity can be non-null on both sides |
| Single grid or multiple grids for non-variants? | 5 | Best answered against real fusion rows, so genuinely wait for Phase 4 |

## Deferred, deliberately

- **Matching stored seqauto text to specimens and extractions** — Phase 2 keeps the client's references
  as text whenever the objects are absent, so the reconciliation pass can land whenever the tracking
  system integration is ready, against real unmatched rows rather than hypothetical ones.
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
