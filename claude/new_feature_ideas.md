# New Feature Ideas for VariantGrid

**Written:** 2026-03-26
**Context:** Based on reading all open issues across SACGF/variantgrid, plus deep codebase exploration of all major apps.

**Project goals framing:** VariantGrid exists to help medical scientists (a) solve disease cases from VCF data, and (b) build a classified variant database that can be mined for deeper insights. Everything below is evaluated against those two goals.

---

## What's already well-served

The platform has impressive breadth — 12+ annotation sources, a DAG-based analysis pipeline, full ACMG evidence tracking, multi-lab discordance detection, ClinVar export, and multi-instance syncing. Features that might look absent but already exist:

- **Internal cohort allele frequency** — `GlobalVariantZygosityCount` already tracks how many samples carry each variant
- **Return-of-results notifications** — already handled via email
- **Compound heterozygote detection** — already implemented in `TrioNode`
- **Analysis templates** — `AnalysisTemplate` class already exists

---

## Feature Ideas

---

### 1. Phenotype-Driven Variant Scoring (Beyond Binary Gene Filtering)

**What:** When a patient has HPO terms, score and rank all candidate variants by how well they fit the phenotype, instead of just filtering to HPO-matched genes. Surface a "phenotype relevance score" column in analysis grids.

**Why this matters:** The current `PhenotypeNode` is binary — it filters to genes associated with the patient's HPO terms. But clinically, scientists want to see all candidates ranked, not just a yes/no gate. A variant in a gene with a strong HPO match should bubble to the top; a variant in a gene with a weak or tangential match should stay visible but de-prioritised. This is how tools like Exomiser work and is one of the biggest productivity improvements for rare disease workup.

**What already exists:**
- HPO terms are extracted from patient phenotype text and linked to ontology terms
- OntologyTerm relationships link HPO terms → genes via OMIM/GenCC
- GenDiseaseClassification has confidence levels (DEFINITIVE down to REFUTED)
- `PhenotypeNode` already queries HPO→gene relationships
- Analysis grids support extra columns from node data

**What needs building:**
- A scoring function: HPO specificity × gene-disease confidence × number of matching HPO terms × variant pathogenicity score
- An optional "phenotype score" column on `SampleNode`/`CohortNode` outputs (not a filter, just a ranking column)
- UI toggle: "rank by phenotype match" reorders the grid
- Patient HPO must be linked to the analysis (currently analysis has no direct patient link — only via Sample → Specimen → Patient)

**Effort:** Medium. The ontology graph is already there. The main work is the scoring formula and wiring patient HPO terms into the analysis context. Exomiser's algorithm is published and could be adapted.

**Benefit:** Very high. Saves hours of manual shortlisting per case. Particularly powerful for trios/families where MOI + phenotype score together give strong signal.

---

### 2. Variant Reclassification Analytics Dashboard

**What:** A set of charts and statistics mining the classification history across all labs, covering: VUS reclassification rates over time, inter-lab agreement rates per gene, time-to-reclassification distributions, most disputed variants/genes, and evidence key usage patterns.

**Why this matters:** The "build a database to mine for insights" goal is currently aspirational. VariantGrid stores rich classification history but provides almost no analytical view of it. A research director or clinician can't answer "which genes in our panel have the most VUS?", "how often do we reclassify?", or "which ACMG criteria does our lab apply most differently from others?" without writing SQL queries.

**What already exists:**
- `ClassificationModification` records every change with timestamp and clinical significance
- Evidence JSON is stored per modification — key usage is queryable
- `DiscordanceReport` tracks inter-lab disagreements
- Lab/organisation hierarchy is in place
- The classification app already has some graphs (graphs page referenced in issues)

**What needs building:**
- A dedicated analytics page (or section on the existing graphs page)
- Pre-computed summary statistics (materialised via Celery beat — expensive to compute on the fly)
- Charts: VUS reclassification rate by month/gene, lab agreement heatmap by gene, evidence key usage bar chart (which criteria are used most/least)
- An "insights" panel on the gene page: "15 labs have classified variants in BRCA2 — 40% pathogenic, 35% VUS, 25% benign"
- Filtering by date range, lab, gene panel, condition

**Effort:** Medium. Most data is in the database. The main work is aggregation queries (with appropriate privacy/permission scoping) and the visualisation layer.

**Benefit:** High for research and quality improvement. Shows the platform is generating knowledge, not just storing it. Valuable for grant reporting, lab audits, and identifying genes that need more consensus work.

---

### 3. Per-Sample QC Dashboard

**What:** A unified quality control summary page per sample/VCF showing: callable site coverage, sex concordance (predicted vs. stated), Ti/Tv ratio, het/hom ratio, and variant count distributions — surfaced alongside a pass/warn/fail status.

**Why this matters:** A negative result from an analysis is only meaningful if the sample passed QC. Currently the `SampleStats` model already calculates chrX het/hom ratio (used for sex inference) and variant counts, but this data is buried and not acted on. A contaminated sample or a sample with wrong genome build can go through the entire analysis pipeline and produce misleading results.

**What already exists:**
- `SampleStats` tracks: SNP/indel counts, het/hom/ref counts, chrX counts (used for `chrx_sex_guess`)
- `SampleStatsPassingFilter` has the same for PASS-filtered variants
- `VCFLengthStats` tracks variant length distributions
- `SeqAutoRun` can import flagstats, FastQC, and coverage QC data from the sequencing backend
- Sex is recorded on `Patient` and `Sample` — mismatch is detectable

**What needs building:**
- A `SampleQCSummary` model or view aggregating all QC signals into pass/warn/fail status
- QC thresholds configurable per `GenomeBuild` and `EnrichmentKit`
- A QC status indicator on the sample page and on the analysis node when that sample is used
- A warning banner in analysis if any source sample has failed QC
- Ti/Tv ratio calculation (not currently computed)
- Integration of seqauto QC data with snpdb sample model (currently separate apps with no bridge)

**Effort:** Medium. The raw data is largely there. The work is standardising thresholds, building the summary model, and surfacing warnings in the analysis UI.

**Benefit:** High. Prevents scientists from acting on bad data. Essential for accredited clinical labs. Particularly important as somatic/tumour samples are added (contamination and purity QC is more critical there).

---

### 4. CandidateSearch: Case-Level Solved/Unsolved Tracking

**What:** Extend `CandidateSearch` with a patient case concept — a named grouping of a patient, their samples, and related analyses — that can be marked solved (with the causal allele recorded) or unsolved. Unsolved cases surface automatically in `REANALYSIS_NEW_ANNOTATION` and `CROSS_SAMPLE_CLASSIFICATION` searches.

**Why this matters:** `CandidateSearch` already finds candidates for reanalysis, but there's no concept of a *case* being solved or unsolved. Every analysis looks the same. A scientist has no way to say "this patient is solved — stop surfacing their variants in reanalysis searches" or to track the lab's overall diagnostic yield.

**What already exists:**
- `CandidateSearch` with `REANALYSIS_NEW_ANNOTATION` and `CROSS_SAMPLE_CLASSIFICATION` types already do the hard work of finding candidates
- `Candidate` model has OPEN/RESOLVED/HIDDEN/HIGHLIGHTED status per variant
- `Sample` → `Patient` linkage exists
- `Candidate.reviewer` and `reviewer_comment` fields are modelled but not exposed in the UI

**What needs building:**
- A `ClinicalCase` model linking a patient, their samples, and analyses — with `CaseStatus` (IN_PROGRESS, SOLVED, CLOSED)
- When SOLVED: record the causal allele(s), inheritance mode, condition
- Exclude SOLVED cases from reanalysis candidate searches (or show them separately)
- Expose the `reviewer` / `reviewer_comment` fields on the candidate grid (currently modelled but commented out)
- A "Cases" list view per lab with filtering by status, gene, condition
- Statistics: solved fraction per lab per year (feeds the analytics dashboard in idea #2)

**Effort:** Medium. `CandidateSearch` infrastructure is there; this adds the case wrapper on top of it.

**Benefit:** High. Gives labs a diagnostic yield metric, prevents re-reviewing solved cases, and turns CandidateSearch from a search tool into a case management workflow.

---

### 5. CandidateSearch: Automatic Triggering and Notifications

**What:** Currently `CandidateSearch` requires a scientist to manually create a new search run. Make `REANALYSIS_NEW_ANNOTATION` and `CLASSIFICATION_EVIDENCE_UPDATE` runs trigger automatically when a new annotation version is imported, and send a notification email summarising what changed.

**Why this matters:** The machinery to find annotation-driven candidates already exists, but it's entirely pull-based — scientists have to remember to go looking. Automatically running searches on annotation import and emailing "5 variants in your active analyses have updated ClinVar significance" transforms it from a tool into an alert system.

**What already exists:**
- `CandidateSearchType.REANALYSIS_NEW_ANNOTATION` — finds analyses with P/LP ClinVar variants under older annotation
- `CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE` — finds classifications with new gnomAD/pathogenicity data
- `CandidateSearchRun` tracks run status and config
- Email notification infrastructure exists
- Annotation version import already triggers other downstream tasks

**What needs building:**
- A signal/hook on `VariantAnnotationVersion` activation that auto-creates `CandidateSearchRun` records for each registered analysis/classification set
- Per-user or per-lab subscription preferences: "notify me when new candidates appear for my analyses"
- A digest email: "New reanalysis candidates found after annotation update vX.Y — 3 analyses affected"
- A "what changed" annotation diff view on each `Candidate` (old ClinVar significance vs. new)
- Scheduled auto-runs for `CROSS_SAMPLE_CLASSIFICATION` (triggered when new classifications are published, not annotation updates)

**Effort:** Medium. The search logic and email infrastructure are in place; the main work is the auto-trigger wiring and the per-user subscription model.

**Benefit:** High. Converts existing functionality from a tool scientists remember to use into a system that proactively surfaces important changes.

---

### 6. Gene Page: Classification Summary Counts

**What:** Add pathogenicity summary counts (total B / LB / VUS / LP / P) above the existing classification groupings datatable on the gene page.

**Why this matters:** The full classification grid already exists on the gene page — labs can see every classified variant. But there's no at-a-glance answer to "how many pathogenic variants does our database have for BRCA2?" A scientist opening the gene page to classify a new variant wants that number immediately, not after scrolling a potentially long table.

**What already exists:**
- Full `ClassificationGrouping` datatable on the gene page (rendered via `{% classification_groupings %}` template tag)
- `ClassificationClassificationBucket` enum with BENIGN / VUS / PATHOGENIC buckets
- `GeneSymbolViewInfo` cached properties pattern — a natural place to add count queries
- The data is all queryable via `GeneSymbol` → `Transcript` → `Allele` → `Classification`

**What needs building:**
- A count query in `GeneSymbolViewInfo` grouped by clinical significance bucket, respecting share-level permissions
- A small summary row above the datatable: e.g. "Pathogenic: 12 · Likely Pathogenic: 8 · VUS: 34 · Likely Benign: 5 · Benign: 11"
- Optionally split by allele origin (germline vs. somatic) since the filter already exists on that page

**Effort:** Small. The query is straightforward; the rendering is a few lines of template. The permission scoping is the most careful part.

**Benefit:** Medium-High. High-value information surfaced with minimal effort. Immediately useful when curating a new variant in a gene the lab has history with.

---

### 7. Classification Completeness Scorer

**What:** For each classification record, compute a "completeness score" indicating what percentage of applicable ACMG criteria have been evaluated (either filled in or explicitly marked N/A). Show this as a progress indicator on the classification form and in list views.

**Why this matters:** Many VUS classifications are VUS not because evidence is lacking, but because curators ran out of time and left criteria blank. A completeness score makes the gap visible, helps lab managers audit quality, and motivates curators to complete their work. It also provides a fairer basis for comparing classifications across labs in discordance reports.

**What already exists:**
- `EvidenceKey` has a well-defined schema with criteria types and applicability
- `EvidenceKeyOverrides` per lab can hide irrelevant keys
- Classification modification stores evidence as JSON — blank fields are detectable
- The classification form already groups evidence by category

**What needs building:**
- A completeness scoring function: (filled criteria + explicitly N/A'd criteria) / total applicable criteria
- Applicability rules: some criteria are only relevant for certain variant types (e.g. PVS1 only for LoF, PS3 only if functional assays exist)
- A progress indicator on the classification form
- A sortable "completeness" column on classification list views
- An optional "completeness threshold" for publication (warn if < X% complete before sharing)

**Effort:** Small-Medium. The evidence key schema and variant annotation data needed for applicability rules are in place. The main complexity is defining sensible applicability rules.

**Benefit:** Medium. Improves database quality over time and gives lab managers visibility into classification backlog. Low risk to implement.

---

## Summary Table

| # | Feature | Effort | Benefit | Serves |
|---|---------|--------|---------|--------|
| 1 | Phenotype-Driven Variant Scoring | Medium | Very High | Case solving |
| 2 | Reclassification Analytics Dashboard | Medium | High | DB mining |
| 3 | Per-Sample QC Dashboard | Medium | High | Data quality |
| 4 | CandidateSearch: Case-Level Solved/Unsolved | Medium | High | Case solving + DB mining |
| 5 | CandidateSearch: Auto-trigger + Notifications | Medium | High | Case solving |
| 6 | Gene Page: Classification Summary Counts | Small | Medium-High | Classification quality |
| 7 | Classification Completeness Scorer | Small-Medium | Medium | Data quality |

---

## Recommended Starting Points

Prioritised by impact-per-effort:

1. **Gene Page: Classification Summary Counts** (#6) — very small effort, immediately useful, no new data model
2. **Classification Completeness Scorer** (#8) — small effort, improves data quality passively
3. **CandidateSearch: Auto-trigger + Notifications** (#5) — medium effort, big improvement over existing functionality for free
4. **Phenotype-Driven Variant Scoring** (#1) — medium effort, very high impact for rare disease workup

Ideas #4 and #5 (both CandidateSearch enhancements) are designed to be built together — the case-level model makes notifications more meaningful, and notifications make the case-tracking workflow feel complete.
