# TSO 500 — order of work

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-07 (revised; earlier phases were planned by Claude Opus 4)
Status: in progress

Sequencing plan for [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431)
and the issues under it. Covers what to build in what order and why, not how to build each piece —
each issue carries its own design.

The test data and its gotchas are in
[`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md).
[`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md) is Phase 8's
reuse design and the one sub-plan still live;
[`1747_specimen_tissue_ontology_plan.md`](1747_specimen_tissue_ontology_plan.md) is the tissue work split
out of #1706. Every other sub-plan was deleted as it landed; what outlived them is in the Done and
Still-open entries below, and the PRs are #1705, #1709, #1712, #1715, #1716, #1718, #1719, #1760,
#1761, #1814 and #1834. The design reasoning that outlived Phase 5 lives in code, in
[`snpdb/gene_level_variants.py`](../../snpdb/gene_level_variants.py).

---

## Done

- **Phase 0 — pipeline fixes** (PR #1709, branch `issue_1506_tso500_pipeline_fixes`).
  `vcf_header_filter_ids()` in `library/genomics/vcf_utils.py` reads declared FILTER IDs off the raw
  header lines, so the DragenExonCNV VCF no longer dies at the first pipe stage on PyVCF's FILTER regex.
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
  `SVTYPE` marks a caller describing an event rather than a segmentation interval, so the DragenExonCNV
  VCF's per-gene `Undetermined` call survives. Where separate VCF rows share a locus and have their depths
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

  `patients/migrations/0016_specimen_tissue_status.py` maps `S` → Affected and `G` → Unknown, `G` being
  indistinguishable from an untouched record, and registers a `ManualOperation` only where a
  `PatientRecord` actually carries a `G`. The patient CSV column is a hard rename to
  `'Specimen Tissue status (Reference/Affected/Unknown)'` — an old spreadsheet fails the import naming the
  missing column, since Germline → Unknown is a loss of meaning rather than a translation.

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

- **Phase 5 — gene fusions as variants and the AllFusions.csv parser (#1506, phase 1)** (PR #1719,
  branch `issue_1506_tso500_phase5_gene_fusions`). The design reasoning lives in
  `snpdb/gene_level_variants.py`'s module docstring — in code rather than here, because a `Variant` with
  no coordinate is a special case every reader of `Variant` eventually trips over.
  Gene-level events get a real `Variant`, anchored on one gene-level contig shared by every build, with
  the gene pair encoded in a symbolic alt and a `GeneFusion` companion model. A `Variant` rather than the
  bare `Allele` the issue originally proposed, because `VariantGeneOverlap` and `CohortGenotype` are both
  `Variant`-keyed — without one, fusions cannot reach gene lists, compound-het detection or any analysis
  node. VEP never sees them; a third annotation pipeline type resolves gene → symbol → release genes
  locally and writes the overlap rows for **both** partners.

  Scoped to ingestion, storage and identity — fusion classification equivalence stays deferred, which is
  what the user-group's "research level project" feedback on #1506 was actually about. Somatic
  classifications already land on `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED`, so nothing regresses.

  Two things the issue's original design assumed otherwise. Every row carries both breakpoints
  (`chr8:128806980` form), so fusions are *not* coordinate-free — what is deferred is the breakend
  representation, not the data, which is captured from day one. And the breakpoints are per-observation:
  `ENTPD3-RPL14` appears three times from one caller with three different 5′ breakpoints, so identity is
  the gene pair and the coordinates live in `CohortGenotype.info` with the rest of the per-row data.

  The one design change that came out of building it: partners are identified by a **`FusionGeneId`** row
  whose pk is the HGNC ID where there is one, and a locally allocated id above 1,000,000 where there is
  not. Clone-based identifiers are routine fusion partners (`RP11-458D21.5`, `AC016683.6` are both in
  the test file), and HGNC-only identity left them unrepresentable. The alt namespaces the two apart
  (`<FUSION:HGNC:nnn>` vs `<FUSION:GENE:nnn>`) because only the first kind means the same thing on
  another deployment — anything leaving the system sends `GeneFusion.canonical_str`, not the number.

  Two changes came out of review on #1719, both worth carrying forward. **One way into the database**:
  the loader writes a genotype VCF (`FORMAT/GT` plus the caller's rows as `INFO`) and hands it to the
  ordinary import path, so `ImportCreateVCFModelForGenotypeVCFTask` builds the VCF/Sample/Cohort and
  `ProcessGenotypeVCFDataTask` COPYs the `CohortGenotype` rows — only the bcftools stages are skipped
  (`upload/vcf/gene_level_vcf_preprocess.py`), since `norm --check-ref=s` reads a base from the fasta a
  gene-level record has no coordinate for. A classification naming `BCR::ABL1` resolves to a variant
  coordinate and goes through the classification upload pipeline the same way. And the **nomenclature is
  VICC's `GENE1::GENE2`**, never the single-hyphen form, which means a read-through transcript; that
  string is written into `VariantAnnotation.hgvs_c` and `hgvs_g`, so the grid columns that read those
  directly show `BCR::ABL1` rather than a blank where g.HGVS is otherwise never blank. Anything parsing
  the vendor's file says so in its name (`DragenTSO500AllFusions…`), while `GeneFusion`, `FusionGeneId`
  and `snpdb.gene_level_variants` stay generic — a classification naming a fusion is a second,
  non-TSO500 source.

  Since a `Sample` belongs to exactly one `VCF`, the file creates its own VCF and Sample rather than
  joining the RNA arm's SpliceVariants sample — Phase 4 ties that sample to the extraction, and Phase
  7's specimen level shows both arms together.

  Follow-ups after the PR, on master directly: `VCFSourceSettings.genome_build`
  (`snpdb/models/models_vcf.py:VCFSourceSettings`) supplies the build for a source whose file has no
  contigs to detect from, so a multi-build deployment configures the fusion CSV's build once per source
  rather than declaring it at every upload; fusions got a variant class in annotation and evidence-key
  options (`annotation/gene_level_annotation.py`); the c.HGVS tag renders the fusion name on the variant
  page; and the AllFusions loader shows in the upload page's file list by default.

- **Phase 7 — analysis grouping node, all levels** (private#223). Extraction level first (PR #1718,
  branch `private_issue_223_tso500_phase7_analysis_grouping_node`), then specimen and patient plus a
  new editor (PR #1814, branch `sample_node_source_levels_223`).
  One entry point gathers every VCF sample belonging to an extraction, specimen or patient, so the DNA
  arm's small-variant, CNV and exon-CNV calls land in one analysis without anyone concatenating files
  outside VG, and the RNA arm joins them at specimen level.

  `SampleNode` gained a `source_level` and extraction/specimen/patient FKs rather than sibling node
  classes; `SampleSourceLevel` (`patients/models_enums.py:SampleSourceLevel`) is the hierarchy's enum
  and `get_sample_group` (`patients/sample_grouping.py:get_sample_group`) is the one resolver for all
  four levels. Patient is the union of `Sample.patient` and `extraction__specimen__patient`, because the
  VCF import carries the extraction down without setting patient while the patient CSV does the
  opposite. Samples resolve **at query time** and each one's filters are ORed as a `pk IN (subquery)` —
  `MergeNode`'s fallback mechanism, factored out of it into `queryset_to_pk_in_q` alongside
  `annotate_and_filter_queryset` so there is one implementation. Each subquery carries only its own
  VCF's genotype join, so `INFO/CN`, `SEGID`, `FORMAT/SM` and Phase 2's preserved FILTER values
  survive — which is what ruled the pre-built cross-VCF cohort out. A single-sample group short-circuits
  to the alias path and produces byte-for-byte the query a sample-level node produces today.

  Per-sample thresholds are a `(node, sample)` child table — sapath#301's different `min_ad` per caller —
  with the node's own fields as the value for any sample without a row, so a sample
  `reconcile_pending_extractions` attaches later inherits rather than gets nothing. VCF FILTER: `PASS`
  is one node-level choice and every other code belongs to its own VCF, since `LowDepth` in a DRAGEN
  small-variant VCF is not the same call as in its CNV VCF —
  `NodeVCFFilter` (`analysis/models/nodes/analysis_node.py:NodeVCFFilter`) stores
  `{"pass": bool, "by_vcf": {...}}` rather than translating filter names into every VCF's codes. The
  `arg_q_dict` cache key folds in every sample's genotype collection and cohort version, and deleting a
  sample bumps nodes grouping on its extraction. Counts are live at group level — there is no single
  cohort's stats row to read.

  Exclusions are reported rather than silently applied: a sample the group leaves out (archived VCF,
  another genome build, no permission) surfaces on the node and through
  `GET /patients/extraction/<pk>/samples?genome_build=`, keyed on the source object rather than the node
  so it answers before the node is saved. Dropping the cohort dropped its `genome_build` FK, so the
  single-build constraint that originally ruled Patient out is gone — an analysis has exactly one build,
  and the node reports what that excluded the same way it reports an archived VCF.

  The editor is one grouped select2 that finds a patient, specimen, extraction or sample, with a tree
  underneath showing the whole patient and the picked row's subtree flagged; the radio on a row is the
  level control, rows outside the pick stay visible but dimmed, and a sample row expands to its VCF's
  own FILTER codes and threshold overrides. On the card the pedigree badge stays at every level (it has
  always drawn the patient) and what the node *is* moves to the class strip via `get_node_strip_label()`;
  chips nest specimen ⊃ extraction ⊃ `VCF ×N`, with an amber chip counting what the group left out.
  `ANALYSIS_SAMPLE_NODE_LEVELS` trims the levels a deployment offers.

  `analysis_templates_tag` (`analysis/templatetags/related_analyses_tags.py:analysis_templates_tag`)
  accepts `extraction=` / `specimen=` / `patient=`, so those pages launch templates — a group counts as
  archived only once every sample it reaches is. `AnalysisTemplate.new_version()` accepts them as
  starting variables, keyed on the level's own FK even though the picker field is called `source`.

  **Per-level add-node menu entries were built and then removed** (`de67cf66c`, 2026-09-04): the menu
  has one Sample node entry and the level is picked in the editor. That is the settled answer to the
  "Add-node menu entries per level" item the extraction-level PR left open.

- **Phase 6, the #1706 / private#2837 half — specimen and extraction grids** (PRs #1761 and #1760,
  branches `issue_1706_specimen_extraction_grids` and `issue_2837_variant_samples_grid_detail`).
  Top-level `specimens` and `extractions` DataTables grids off the patients menu, both through
  `filter_for_user`; the extraction grid's sample count is annotated over `Sample.filter_for_user`
  rather than the raw related set, since seeing an extraction does not mean seeing what was sequenced
  off it. The samples grid links specimen and extraction (the specimen reference alone cannot tell the
  DNA and RNA arms apart). The variant page's samples grid gained a click-to-expand row carrying patient
  code, specimen, tissue, tissue status, collection date, extraction and nucleic acid, each linked,
  with the same three as hidden columns so the search box filters by them — the identifiers join onto
  the existing single query, and detail is only given for patients passing `Patient.filter_for_user`,
  patient code only, never name. Patient falls back to the specimen's patient, since the sequencing
  pipeline links a sample to an extraction without touching `Sample.patient`. The tissue dropdown is
  excluded from the specimen forms until #1747 gives `Tissue` a way to be created.

- **Phase 8, sapath#246 — Classify & Report tab** (PR #1834, branch `246_classify_and_report`), which
  also lands the first half of #444. A tab on the sample page (and the patient page, unioned over the
  patient's samples) that turns outstanding classification tags into classifications, then builds a
  multi-variant report from them. `analysis/classify_report.py` and
  `analysis/views/views_classify_report.py`, in the analysis app because it is built on `VariantTag`.

  `Tag.requires_classification` marks the tags that feed the queue (`RequiresClassification`,
  `SomaticReportable`; the latter also gets the somatic allele-origin bucket). `VariantTag.sample` is the
  **study's proband** from `AnalysisNode.get_proband_sample()` on the node it was tagged in — carrying
  the variant is deliberately not what decides ownership, since a relative HET for the proband's variant
  does not need their own classification. Nullable forever: filled silently where the node knows the
  proband, else picked at classification time. Existing taggings are backfilled by
  `one_off_backfill_variant_tag_sample`, registered as a `ManualOperation`. A to-do tagging is
  **resolved against a classification rather than deleted** (`VariantTag.resolved_classification`), so it
  stays as the record of what was flagged and what it turned into; withdrawing the classification puts
  the to-do back. Resolution is automatic when the classification is of the tagging's own sample or the
  analysis is about one person (`analysis/variant_tag_operations.py:classification_resolves_tag`);
  otherwise the row shows the link and a Clear tag button so the scientist says which person it was for.

  The Classify dialog lists the lab's previous classifications of the allele with condition, curated
  date and summary — which one applies is always the scientist's call, so nothing is highlighted.
  "Apply to this sample" creates the record and copies that record's consensus; "New classification"
  opens the full form. `ClassificationReport`
  (`classification/views/classification_export_report.py:ClassificationReport`) accepts a list of
  classifications and an explicit template, with `classifications`, `gene_groups` and a case header in
  the context alongside the unchanged single `record`, so existing single-variant templates keep working.

## Where things stand

The five test files from `ed5e15a33` — one de-identified run, one specimen, DNA and RNA arms, under
`upload/test_data/tso500/ExampleSample_2600000001/`:

| File | State |
|---|---|
| hard-filtered.vcf | loads, 93 records |
| cnv.vcf | loads, 16 records — 9 copy-neutral rows skipped and counted; `SM` surfaced by Phase 6 |
| SpliceVariants.vcf | loads, 17 records — VAF derived from `ALTDEDUP`/`REFDEDUP`, `LowUniqueAlignments` preserved |
| DragenExonCNV.vcf | loads, 2 records, against a `genome_build` declared at upload or on its `VCFSourceSettings` row |
| AllFusions.csv | loads, 33 rows -> 31 fusion variants (2 gene pairs seen more than once); build from `VCFSourceSettings.genome_build` |

Once Phase 4 has tied their samples to `2600000001C` and `2600000001B`, the DNA arm's three files reach
one analysis through an extraction node and all five through a specimen node, each row carrying its own
VCF's values. None of that has yet been run end to end on a real deployment — see the manual pass
below.

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
      │              │              P7  private#223 extraction ✔ specimen ✔ patient ✔
      │              │
      │         P8  sapath#246 ✔ ──► copy_consensus audit ──► #1419 ──► #444 remainder
      │
  P5  #1506 GeneFusion ✔ ──► AllFusions parser ✔
      │
  P6  #1706 grids ✔  private#2837 ✔  ──► #1558 kind badge ✔ copy number ✔ fusion columns + filter ✔
                                          #1747 tissue as UBERON
```

Nothing left has a hard dependency on anything else; what remains is ordered by value.

---

## What remains

1. **Phase 8 remainder** — the `copy_consensus` audit, #1419's gene-level copy, and #444's report
   template and specimen-page launch. Below.
2. **The manual pass over a real GRCh37 deployment** — owed by Phases 2, 4, 5 and 7 and never run.
   One session, listed under each phase below.
3. **Client work** — `variantgrid_api#20`, deferred to batch with the other client issues.
4. **Chase TAU** for the missing files. Long lead time; below.
5. **#1717** — deployment-wide threshold defaults per VCF source and panel. Independent.
6. **#1747** — tissue as a UBERON term. Independent; PR #1761 hid the tissue dropdown until it lands.

## Still open from Phase 0 — chase TAU for the missing files

Long lead time, so worth pushing now:

- `Logs_Intermediates/Gis/<sample>/<sample>.abcn_annotated.vcf` + `<sample>.abcn_genes.tsv` — per-gene
  absolute and minor copy number, the latter being where gene-level LOH comes from. The only genuinely
  CVO-only data, and a VCF, so it is in ingestion scope. Illumina names the path but publishes no
  INFO/FORMAT spec, so it cannot be mocked up.
- MetricsOutput.tsv / the run-level TMB summary — where the numbers Phase 4's `SpecimenMeasure` API
  takes are transcribed from, so it is now the client's need rather than VG's.
- A run with a real BRCA1/2 large rearrangement. Both DragenExonCNV records in the test data are
  constructed, and no published TSO 500 output has a populated one either. Until one arrives the field
  spelling in that file stays provisional.

## Still open from Phase 1

The patient CSV columns are all `SPECIMEN_*`, so a CSV can name one extraction per specimen but not the
DNA and RNA arms of one block. Worth revisiting when a consumer needs it.

## Still open from Phase 2

The manual checklist on #1711 wants running against a real GRCh37 deployment: PR #1712 measured the
record counts and VAFs through `write_cleaned_vcf_header` → `vcf_clean_and_filter` → `bcftools norm`,
not through a full import with annotation.

The `source` strings a client sends (`DRAGEN TSO500 SmallVariant`, `DRAGEN TSO500 CNV`) are documented
in [`upload/test_data/tso500/README.md`](../../upload/test_data/tso500/README.md) but nothing is keyed
on them yet — the SpliceGirl mapping comes off the header and the copy-neutral skip is a general rule.
They become part of VG's configuration contract the moment a `VCFSourceSettings.source_regex` matches
one, so they want to stay stable from the first client.

`SEGID=MYCL1` resolving to `MYCL` is unconfirmed against a real database. NCBI carries the alias, and
`GeneSymbolMatcher.get_gene_symbol_id_and_alias_id` (`genes/gene_matching.py:45`) is the resolver, but
Phase 6's gene-symbol item rests on that assumption holding. Phase 5's fusion parser goes through the
same resolver for `SEPT14` → `SEPTIN14`, so one check against a real database covers both — and if the
alias is missing, a fusion partner still imports, just under a local `GENE:` id rather than its HGNC one.
[`1669_release_gene_matcher_alias_chaining_plan.md`](1669_release_gene_matcher_alias_chaining_plan.md)
changes alias resolution to single-hop, so do the check after that lands.

**Clients send a build's own name (`GRCh37`), not an alias.** `GenomeBuild.get_name_or_alias("hg19")`
raises `MultipleObjectsReturned` rather than `DoesNotExist`, so a declared build that will not resolve
is a 400 rather than a guess. `hg19` happens to resolve here because that build is disabled, but
`GenomeBuild.enabled` is per-deployment DB state, so the check stays. Documented in the API schema,
`import_vcf --genome-build` and the test-data README; the lab client Phase 4 unblocks should follow it.

## Still open from Phase 3

A deployment whose patient CSV ever populated the old Germline/Somatic column gets a manual task out of
`patients/migrations/0016_specimen_tissue_status.py` asking a human to decide which of those specimens
really were Reference. Nothing registers on deployments that never used the column.

`Tissue` has no way to be created, which is #1747 and its own plan.

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

## Still open from Phase 5

The fusion loader's tests assert the VCF it writes — records, sample column, genotype, merged
observations — rather than driving the pipeline steps by hand; what happens to that VCF afterwards is
the ordinary genotype import path the existing VCF suite covers. So a real end-to-end upload of
AllFusions.csv joins Phases 2 and 4 on the manual pass over a real deployment, where the
`SEPT14` → `SEPTIN14` alias check also lives.

## Still open from Phase 7

One run joins the manual list: an extraction node over the real TSO 500 run returning the DNA arm's
three callers in one grid, each row carrying its own VCF's `INFO`/`FORMAT`, different `min_ad` per
caller, ×3 on the canvas; then a specimen node over the same run pulling the RNA arm's splice and fusion
rows in beside them. PR #1718's and #1814's evidence stops at the unit suites over synthetic VCFs, and
PR #1814 notes the editor's tree and the card chips were not checked in a browser.

Deployment-wide threshold defaults per VCF source and per panel are #1717 rather than part of this
phase — the per-node values built here run either way; #1717 only changes what a new node starts from.

## Phase 6 — #1558

Landed. Designed in [`1558_non_variants_on_grids_plan.md`](1558_non_variants_on_grids_plan.md), which
is the spec and records what each phase touched - the kind badge, the fusion row expansion and export,
copy number per sample, the fusion calls column and the Effect node's structural filter.

Per the issue's own comment table, SV already worked in the grid and CNV worked as SV — what was
missing was surfacing `CohortGenotype.info["CN"]` and TSO 500's `FORMAT/SM` linear copy ratio, which
importer v21+ already kept in the format JSON blob but which was not queryable. That is now
`VCF.copy_number_field` plus a read of the JSON at query time.

cnv.vcf's `SEGID` gene symbol came down from Phase 2 with it, and was deliberately left unsurfaced:
it is whatever the caller wrote — `MYCL1` where the rest of the pipeline says `MYCL` — and the gene and
overlapping-symbol columns already say which gene a segment hits, resolved properly.

Fusions were the genuinely new case. **One grid, not several** — a fusion is a row like any other, so
compound-het detection, gene lists and every downstream node keep working over the same result set
rather than needing a per-kind union. That is also what Phase 5's choice of a real `Variant` over a bare
`Allele` was for. Phase 5 left them grid-ready rather than grid-complete; what this phase added is the
kind badge that says a row is one, the partners and direction in the expanded row, the Fusion calls
column (breakpoints, caller and read counts out of `CohortGenotype.info`, one entry per row the caller
wrote), the export and All Variants paths they used to fall out of, and a structural filter on the
Effect node.

Finding a fusion by a **single** partner — `ROS1` returning every fusion it takes part in — is the gene
page: it already joins `VariantGeneOverlap`, which the gene-level annotation run writes for both
partners, and the badge now makes those rows recognisable among the small variants. Phase 5's search
receiver still resolves only a full pair (`BCR::ABL1`, `CD74-ROS1`, aliases applied) and stays
lookup-only, minting no `GeneFusion` or `FusionGeneId` on whatever a user types.

## Phase 8 — what is left (#1419, then #444)

Designed in [`somatic_curation_reuse_issue_1419_plan.md`](somatic_curation_reuse_issue_1419_plan.md),
which is the spec. Somatic reporting sees the same variants over and over, so the phase is about
reusing prior curation. sapath#246 landed the entry point: the queue, the dialog listing prior
classifications, and "Apply to this sample" copying consensus. Still to do, in order:

- **The `EvidenceKey.copy_consensus` audit.** 142 keys carry it, set for germline ACMG work and never
  reviewed against somatic, and several of them copy one patient's tumour measurements or report
  narrative onto another's. "Apply to this sample" copies consensus today, so this audit is now the
  most urgent item in the phase — until it is done, the tab can carry a prior patient's data across.
- **#1419** — gene-level copy consensus, deliberately as a copy rather than the first-class gene/disease
  object the issue also proposes. The object needs the disease axis settled first; the copy reuses
  machinery that already works and owns exactly the key set the object would. The human picks the
  source record because AMP tiering is gene *and* tumour type and the phenotype data cannot make that
  call — which is what the #246 dialog already does for the same allele.
- **#444 remainder** — `ClassificationReport` takes a list and the context carries `gene_groups`, but the
  vue report template still wants an array-aware version, Phase 4's TMB/MSI `SpecimenMeasure`s are not
  yet pulled into the report, and the tab launches from sample and patient pages but not the specimen
  page. Linked-classification semantics (compound het / ID linkage) were left out of #1834 as
  selection at report time covers the reporting need without new data.

---

## Parallelism

Three independent streams: Phase 6's columns and fusion filter (`snpdb`, `analysis` grids), Phase 8's
audit and #1419 (`classification`), and the manual pass over a real deployment, which needs a person
with a GRCh37 box and the API token rather than a developer. #1747 and #1717 slot in anywhere.

With one person: the `copy_consensus` audit first, since the #246 tab now exercises it on real
patients; then the manual pass, which validates everything landed so far before more is built on it;
then Phase 6; then #1419 and #444.

## Decisions settled in code

Phase 5's identity question — `FusionGeneId` gives every partner an identity, HGNC where there is one
and a local id where there is not, so no row is parked. Phase 6's single-vs-multiple grid question —
one grid, so fusions reach comp-het and every downstream node. Phase 7's — per-sample threshold rows in
the editor, restrict to the analysis build and report the exclusions, `Sample.extraction` left as a
single FK, live counts at group level, and one Sample node menu entry with the level picked in the
editor. Phase 8's — a tagging is resolved against its classification rather than deleted, the tagging's
sample is the study's proband rather than a carrier, and which prior record to copy is always the
scientist's choice.

## Deferred, deliberately

- **Per-gene LOH as a specimen measure** — #1559 lists it, but it is per-gene rather than per-specimen and
  the vendor emits it in a VCF, so it belongs with abcn_annotated.vcf below rather than in Phase 4.
- **Fusion equivalence and discordance** — the research-level part of #1506. Somatic classifications get
  `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED` today, so nothing regresses by waiting.
- **Gene/disease curation as a first-class object** — #1419's schema option, and the only thing that can
  express "reviewed on this date, due for review". Needs the disease axis settled; Phase 8 does the
  reuse as a copy in the meantime, over exactly the key set such an object would own.
- **Breakend representation and BND VCF export** — #1506 phases 2 and 3. The breakpoint values themselves
  arrive in AllFusions.csv and Phase 5 stores them, so this is a read-side change when a consumer appears.
  Whether VEP parses BND ALT syntax at all is unverified — a 20-minute experiment before anyone designs
  around either answer — but what VEP would give is per-position feature overlap rather than frame or
  domain analysis, and the caller already reports that as `Gene A/B Location`.
- **Coordinate-free gene-level CNV** — designed alongside Phase 5 and purely additive to it:
  `<AMP:HGNC:nnn>` / `<LOSS:HGNC:nnn>` on the same gene-level contig, same anchor, same annotation run —
  one enum value and one alt prefix. Only for a caller that reports a gene-level event with no
  coordinates at all (the CombinedVariantOutput "JAK2 amplification (5 copies)" style); cnv.vcf's
  `<DUP>`/`<DEL>` carry real coordinates and stay structural variants. Copy number stays
  observation-level in `CohortGenotype.info` rather than in the alt, because labs use different
  amplification thresholds and per-count identity would stop two labs ever agreeing on "JAK2
  amplification". Build it when a file needs it.
- **abcn_annotated.vcf / gene-level LOH** — cannot start until TAU supplies a file; the format is
  undocumented and not usefully mockable.
- **DragenExonCNV exact field spelling** — provisional until a run with a real large rearrangement
  arrives. Note the CombinedVariantOutput reports these as `<LOSS>` where the VCF header declares `<DEL>`.
- **CombinedVariantOutput.tsv** — settled as not worth ingesting. It is exactly a filtered subset of the
  individual files (`PASS` small variants, `PASS` non-reference CNV) and discards 148 of 149 fusion calls
  and every splice call.
- **Fold-change → DEL/DUP conversion** (sapath#304's open question) — moot for v2.6.2, which emits
  `<DUP>`/`<DEL>` directly with `SM` as the linear copy ratio. The `cnv_tsv_to_vcf.py` command in the
  sapath repo predates that.
