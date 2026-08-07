# Data collection plan — VariantGrid paper (clinical + research usage mining)

_Companion to `data_mining_opportunities.md` (the "why/what figures") and `scripts/` (the "how").
This file specifies the concrete datasets, the model/field each comes from, the governance rules,
and the VG3→VG4 approach. Scripts are read-only and produce AGGREGATE CSVs only._

## 0. Governance & safety (read first — non-negotiable)
- **Approvals:** clinical-usage analysis of SA Pathology data needs the appropriate
  ethics/quality-assurance pathway (HREC approval or a QA/audit exemption per SA Health policy).
  Confirm before any numbers leave the clinical environment. This is the long pole — start now.
- **Aggregate only.** Every output is a count/rate/distribution. **No patient-level rows leave.**
- **No identifiers exported.** Never emit: patient IDs, sample names, usernames, condition
  free-text, HGVS/coordinates tied to a person. Users are reported as *counts of distinct users*,
  never identities.
- **Small-cell suppression.** Any aggregate cell with count < `MIN_CELL` (default 5) is suppressed
  (reported as `<5` or dropped). Prevents re-identification from rare combinations.
- **Read-only.** Scripts only read (counts/aggregations). No writes, no migrations. Prefer running
  against a read replica or a restored backup, not live prod.
- **Claude does not run on the clinical environment.** These scripts are written here for *you* to
  review and run. Treat them as drafts to verify in-context.

## 1. VG3 vs VG4 approach
- Clinical data lives on **`vg3_sapath_prod`** (in production since 2020, until the VG4 upgrade in
  a few months). Research data on variantgrid.com.
- **Recommended: migrate the data up and run the VG4 versions of these scripts.** The scripts here
  target **current (VG4/master) models**, so running post-migration avoids schema-version drift.
- After migration, **sanity-check 2-3 aggregate numbers** against known VG3 values (e.g. total
  samples, total classifications) to confirm the migration preserved what you're counting.
- Where a field may not exist or backfill cleanly in older data (e.g. timestamps on very early
  records), the scripts mark it `# VERIFY` and degrade gracefully.

## 2. Datasets → paper claim → source (each becomes a CSV)

| # | Dataset (CSV) | Backs which claim/figure | Source model · field |
|---|---------------|--------------------------|----------------------|
| 1 | Ingestion over time (VCFs, samples / month) | Growth; "always-on ingestion" | `snpdb.VCF.date`; `Sample` via `vcf.sample_set`; `Cohort.sample_count` |
| 2 | Genotype:variant ratio per build | SQL multi-sample packing scales; managed annotation | `snpdb.Variant` count; `CohortGenotype` rows × `CohortGenotypeCollection.num_samples` |
| 3 | Annotation versions + diffs (added/modified/removed) | Versioned annotation + diff-driven re-analysis (VG-native) | `annotation.VariantAnnotationVersion`; `VariantAnnotationVersionDiff.num_*`; `VersionDiffFromToResult` |
| 4 | ClinVar/clin-sig churn between versions | What actually changes between annotation versions | `VersionDiffFromToResult` (value_from, value_to, count) for the clin-sig column `# VERIFY column name` |
| 5 | Classifications over time (by sig, share level) — LIGHT | "we're substrate beneath Shariant"; reanalysis context | `classification.Classification.created`, `clinical_significance`, `share_level`, `lab` |
| 6 | Significance-change flags / month | Reclassification-over-time happens in practice | `flags.Flag` where `flag_type__id='classification_significance_change'`, `.created` |
| 7 | Analysis usage over time + nodes/analysis | Real interactive use of the node-graph | `analysis.Analysis.created`; node counts via `AnalysisNode` MTI subclasses |
| 8 | Node-type usage histogram | Which node types people actually use (live node-graph value) | each `AnalysisNode` subclass `.objects.count()` |
| 9 | Template runs vs ad-hoc; template modification rate | "exploration vs fixed clinical pipeline" thesis | `AnalysisTemplateRun` count; nodes changed after template run |
| 10 | Variant tagging (by tag type, over time) | Curation behaviour; tag→classify funnel | `analysis.VariantTag.tag_id`, `.created` |
| 11 | Already-classified-variant recurrence | Internal population DB / matching value | alleles with a `Classification` ∩ samples carrying them (`VariantZygosityCount`/`CohortGenotype`) — ADVANCED, verify perf |
| 12 | Search activity (best-effort) | Search/lookup usage | `eventlog.Event` / `ViewEvent.args` — NO dedicated search log; best-effort only |
| 13 | Cross-deployment usage contrast | Configurability (one codebase → 4 deployments) | run 1-10 on each deployment; compare profiles |

## 3. Output format
- One CSV per dataset, written to an `--output-dir`. Header row + aggregate rows only.
- Time series bucketed by **month** (also provide yearly rollups). Columns kept minimal.
- A `_run_manifest.csv` records: script version, date run, deployment name, row counts emitted,
  `MIN_CELL` used, and any collectors that were skipped/failed — so the provenance of every figure
  is auditable.

## 4. Highest-value first (if time-boxed)
Do **1, 2, 3/4, 7/8** first — they back the *architecture* claims (ingestion scale, versioned
annotation + diffs, live node-graph use) that are the paper's spine. Then 6 and 10 (reanalysis +
curation behaviour). 11 (recurrence) is high-value but heaviest — do once basics work. 5 stays
deliberately light (Shariant's turf). 12 is best-effort.

## 5. What I (Claude) can and cannot do here
- **Can:** write/maintain the read-only scripts (`scripts/`), define the aggregations, document how
  to run them, and help interpret the resulting CSVs once you bring the *aggregate* outputs back.
- **Cannot / will not:** run anything against the clinical environment. You run the scripts; you
  review every query before running; only de-identified aggregates come back out.
