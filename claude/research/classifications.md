# Classification App: Detailed Technical Reference

## Purpose

The `classification` app is the largest and most complex app in VariantGrid. It implements a full **variant classification workflow** for clinical and research genomics:

- Curate genetic variant interpretations using ACMG/AMP and custom evidence schemas
- Manage versioned, audited classification records with multi-lab sharing controls
- Detect and resolve conflicting classifications between labs (discordance)
- Submit classifications to ClinVar
- Match variant conditions to standardized disease ontology terms
- Import/export classifications in multiple formats (CSV, REDCap, JSON, ClinVar XML)

In clinical genomics, a "classification" is a formal assessment of whether a genetic variant is: **Benign (B)**, **Likely Benign (LB)**, **Variant of Uncertain Significance (VUS)**, **Likely Pathogenic (LP)**, or **Pathogenic (P)** — or for somatic variants, an AMP Tier (I–IV).

---

## Directory Structure

```
classification/
├── admin/                            # Django admin interfaces (12+ admin classes)
├── autopopulate_evidence_keys/       # Auto-fill evidence from variant annotation
├── enums/                            # Enums (clinical significance, share levels, discordance)
├── management/commands/              # Django management commands
├── migrations/                       # Database migrations
├── models/                           # Core models (30+ files)
├── signals/                          # Django signal handlers
├── tasks/                            # Celery async tasks (4 task files)
├── templates/classification/         # Django HTML templates
├── templatetags/                     # Custom template tags
├── tests/                            # Unit and integration tests
├── utils/                            # Utilities (ClinVar matcher, HGVS, etc.)
├── views/                            # Views and API endpoints (28+ view modules)
├── forms.py                          # Django forms
├── serializers.py                    # DRF serializers
├── urls.py                           # URL routing (100+ patterns)
├── apps.py                           # App config
├── classification_changes.py         # Change tracking for activity log
├── classification_import.py          # Import processing
├── classification_stats.py           # Statistics computation
├── criteria_strengths.py             # ACMG criteria strength handling
├── evidence_key_rename.py            # Evidence key refactoring utilities
└── variant_card.py                   # Summary card generation
```

---

## Core Models

### `Classification` (`models/classification.py`)

The primary record for a variant classification. Key fields:

| Field | Type | Description |
|-------|------|-------------|
| `lab` | FK → Lab | Which lab created this classification |
| `variant` | FK → Variant | The genetic variant being classified |
| `sample` | FK → Sample (nullable) | Optional sample context |
| `genome_build` | FK → GenomeBuild | Target genome build |
| `allele` | FK → Allele (nullable) | Cross-build allele link |
| `submission_source` | CharField | How it was created (API, form, import, etc.) |
| `evidence` | JSONField | All evidence key/value pairs (see Evidence System) |
| `clinical_significance` | CharField | Cached current significance (B/LB/VUS/LP/P) |
| `share_level` | CharField | Visibility (user/lab/institution/public) |
| `withdrawn` | BooleanField | Soft-deleted flag |
| `user` | FK → User | Creator |
| `last_edited_version` | FK → ClassificationModification | Current version |
| `last_published_version` | FK → ClassificationModification | Last published snapshot |

**Key behaviors:**
- All data changes create a new `ClassificationModification` (versioned immutably)
- Publishing broadcasts to a share level; can't be revoked, only withdrawn
- Auto-populates evidence fields from variant annotation on creation
- Linked to `ClinicalContext` for discordance tracking

---

### `ClassificationModification` (`models/classification.py`)

Immutable snapshot of a classification at a point in time:

| Field | Type | Description |
|-------|------|-------------|
| `classification` | FK → Classification | Parent classification |
| `evidence` | JSONField | Full evidence blob at this point |
| `is_last_edited` | BooleanField | Most recent edit |
| `is_last_published` | BooleanField | Most recent published version |
| `publish_level` | CharField | Share level when published |
| `previous` | FK → ClassificationModification | Linked list of history |
| `created` | DateTimeField | Timestamp |
| `user` | FK → User | Who made this change |

Used for diffs, history views, reverting, and regulatory audit trails.

---

### `EvidenceKey` (`models/evidence_key.py`)

Defines the schema for all evidence fields on a classification. Each `EvidenceKey` represents one field in the evidence JSON blob.

| Field | Type | Description |
|-------|------|-------------|
| `key` | CharField | Unique string identifier (e.g., `"clinical_significance"`) |
| `value_type` | CharField | Type of value (see EvidenceKeyValueType enum) |
| `evidence_category` | CharField | Grouping category |
| `mandatory` | BooleanField | Must be populated before sharing |
| `options` | JSONField | Predefined choices for select fields |
| `default_crit_evaluation` | CharField | ACMG strength default |
| `max_share_level` | CharField | Restricts visibility |
| `hide` | BooleanField | Hide from UI |
| `immutable` | BooleanField | Can't be changed after creation |
| `namespaces` | ManyToMany | Namespace groupings |
| `order` | IntegerField | Display ordering |
| `help_text` | TextField | Field description |

**EvidenceKeyValueType** options:
- `free` — Free text entry
- `textarea` — Multi-line text
- `select` — Single choice from options
- `multiselect` — Multiple choices
- `boolean` — Yes/No
- `date` — Date picker
- `age` — Age with units
- `criteria` — ACMG criteria strength (BA1, BS1, etc.)
- `user` — User selector
- `unit` — Float 0.0–1.0
- `integer` — Integer
- `float` — Floating point
- `phenotype` — HPO phenotype term

---

### `EvidenceMixin` (`models/evidence_mixin.py`)

Base class used by `Classification` and `ClassificationModification`. Provides:

- Evidence CRUD (get/set/delete individual keys)
- ACMG criteria evaluation (computing classification from evidence)
- Patch operations (merge evidence dicts)
- Validation logic
- Blob structure with `value`, `note`, `explain`, `db_refs`, `validation`

**Evidence blob structure per key:**
```json
{
  "clinical_significance": {
    "value": "LP",
    "note": "Met PM1, PM2, PP2, PP3",
    "explain": "See literature review",
    "db_refs": [{"db": "PubMed", "id": "12345678", "url": "...", "summary": "..."}],
    "validation": [{"code": "warning", "message": "..."}]
  }
}
```

---

## The Evidence Key System

Evidence keys are the core schema for classifications. They define what data fields exist and how they behave.

### Special Evidence Keys (`enums/SpecialEKeys`)

These keys have hardcoded behavior in the system:

| Key | Description |
|-----|-------------|
| `variant_coordinate` | Genomic coordinate (chr:pos ref>alt) |
| `c_hgvs` | Coding HGVS notation |
| `g_hgvs` | Genomic HGVS notation |
| `p_hgvs` | Protein HGVS notation |
| `clinical_significance` | Main classification (B/LB/VUS/LP/P) |
| `somatic_clinical_significance` | Somatic tier |
| `allele_origin` | Germline / Somatic / Unknown |
| `condition` | Disease/phenotype text |
| `gene_symbol` | Gene symbol |
| `transcript` | Transcript accession |
| `clingen_allele_id` | ClinGen Allele Registry ID |
| `gnomad_af` | gnomAD allele frequency |
| `spliceai` | SpliceAI prediction score |
| `mode_of_inheritance` | Inheritance pattern |
| `zygosity` | Observed zygosity |
| `sample_id` | Sample identifier |
| `patient_id` | Patient identifier |
| `age`, `sex` | Patient demographics |
| ACMG criteria | `pvs1`, `ps1`–`ps4`, `pm1`–`pm6`, `pp1`–`pp5`, `ba1`, `bs1`–`bs4`, `bp1`–`bp7` |

### Lab Configuration of Evidence Keys

Labs can customize evidence keys via `EvidenceKeyOverrides`:
- Override visibility, options, defaults
- Enable/disable namespaces (e.g., disable ACMG criteria for somatic model)
- Configuration merges: Organization → Institution → Lab level

### `EvidenceKeyMap`

Singleton cache of all evidence keys with lab/org overrides applied. Used throughout the app for key lookup, option resolution, and namespace filtering.

---

## Share Levels

Classifications have a `share_level` controlling visibility:

| Level | Description |
|-------|-------------|
| `user` | Visible only to the creating user |
| `lab` | Visible to all lab members |
| `institution` | Visible to the institution |
| `logged_in_users` | Visible to all authenticated users |
| `public` | Publicly visible (no login required) |

Fields themselves can also have `max_share_level` to restrict sensitive data (e.g., patient info) even when the classification is public.

---

## Versioning System

### Publish Flow

```
[Draft evidence edits] → patch → [ClassificationModification (is_last_edited)]
                ↓
         [Publish at share level]
                ↓
         [ClassificationModification (is_last_published)] ← locked snapshot
```

- Publishing is one-way (can't revert to draft)
- Multiple draft edits before publishing accumulate in history
- Withdrawing soft-deletes the whole classification (visible to lab, hidden to others)
- Full history accessible via linked list of `ClassificationModification` records
- Diffs generated between any two versions

---

## Clinical Context and Discordance

### `ClinicalContext` (`models/clinical_context_models.py`)

Groups all classifications for the **same allele** in order to detect discordance. Each allele can have multiple clinical contexts (e.g., one per disease condition).

| Field | Type | Description |
|-------|------|-------------|
| `allele` | FK → Allele | The variant being assessed |
| `name` | CharField | Context name (often the condition) |
| `condition` | FK → ConditionTextMatch | Linked ontology term |
| `pending_cause` | CharField | Why recalculation is pending |
| `discordance_status` | CharField | Current discordance level |

When a classification is published or changed, `ClinicalContext.recalculate()` runs:
- If all classifications agree → Concordant
- If classifications disagree → Discordant → triggers `DiscordanceReport`

### `DiscordanceReport` (`models/discordance_models.py`)

Tracks a detected discordance between labs:

| Field | Type | Description |
|-------|------|-------------|
| `clinical_context` | FK → ClinicalContext | Context where discordance was found |
| `resolution` | CharField | Ongoing / Concordant / Continued Discordance |
| `cause_text` | TextField | Admin explanation of resolution |
| `closed_by` | FK → User | Who resolved it |
| `closed` | DateTimeField | When resolved |
| `report_started_date` | DateField | When discordance was detected |

### `DiscordanceReportClassification` (`models/discordance_models.py`)

Links specific classifications to a discordance report, with their significance at time of detection.

### `DiscordanceReportTriage` (`models/discordance_models.py`)

Per-lab response to a discordance report:
- Labs can acknowledge, note planned review, or mark as acceptable
- Tracks which lab has triaged and what action they plan

### Discordance Levels (Enum)

```
No entries → Single submission → Concordant → Discordant
```

Discordant sub-types exist for B vs P extremes vs VUS ranges.

---

## Condition Matching System

### `ConditionText` (`models/condition_text_matching.py`)

A normalized free-text condition string entered by a lab. Unique per lab + text. Tracks how many classifications reference it.

### `ConditionTextMatch` (`models/condition_text_matching.py`)

Hierarchical mapping of condition text to ontology terms:

```
ConditionText (root)
├── ConditionTextMatch (gene-level) — matches for this condition + gene
│   └── ConditionTextMatch (classification-level) — specific classification override
```

Each match can reference:
- One or more ontology terms (MONDO, OMIM, HPO, Orphanet)
- Mode of inheritance refinement
- Multi-condition logic: `AND` / `OR` / `NOT_DECIDED`

### Matching Algorithm (`utils/ontology_matching.py`)

1. Normalize text (lowercase, strip punctuation)
2. Search ontology via fuzzy matching
3. Walk ontology relationships (ancestors/descendants)
4. Support multi-condition strings (comma-separated, "Condition A and Condition B")
5. Manual curation workflows for unmatched terms

---

## ClinVar Integration

### `ClinVarAllele` (`models/clinvar_export_models.py`)

Wraps an `Allele` for ClinVar submission purposes. Tracks all ClinVar exports for that allele.

### `ClinVarExport` (`models/clinvar_export_models.py`)

Represents one unique ClinVar submission: a combination of Lab + Allele + Condition. Key fields:

| Field | Type | Description |
|-------|------|-------------|
| `clinvar_allele` | FK → ClinVarAllele | The allele being submitted |
| `condition` | FK → ConditionTextMatch | Resolved condition term |
| `classification_based_on` | FK → ClassificationModification | Source classification version |
| `status` | CharField | Pending / In Batch / Submitted / Error |
| `scv` | CharField | ClinVar SCV accession returned after submission |

### `ClinVarExportBatch` (`models/clinvar_export_models.py`)

Groups multiple `ClinVarExport` records for a single batch API submission.

### `ClinVarExportSubmission` (`models/clinvar_export_models.py`)

Immutable record of data sent and response received for one submission attempt.

### ClinVar Export Convertor (`utils/clinvar_export_convertor.py`)

Maps VariantGrid evidence keys to ClinVar's XML/JSON submission format:
- Maps clinical significance values
- Converts condition terms to MedGen IDs
- Handles inheritance mode mapping
- Generates assertion method references
- Packages citation/pubmed references

### ClinVar Matching (`utils/clinvar_matcher.py`)

For labs joining Shariant that had prior ClinVar submissions: matches existing ClinVar records to newly created Shariant classifications based on variant, condition, and lab.

---

## Import System

### Upload Pipeline

```
User uploads file (Excel/CSV/ClinVar XML/etc.)
        ↓
UploadedClassificationsUnmapped — stores file, status = Pending
        ↓
[Celery Task] ClassificationImportMapInsertTask
        ↓
omni_importer — external tool maps file to VG JSON format
        ↓
BulkClassificationInserter — creates/updates Classification records
        ↓
VariantResolver — matches HGVS/coords to database Variant records
        ↓
[VCF pipeline] — new variants inserted if not found
        ↓
[Liftover] — coordinates lifted over to other genome builds
```

### `ImportedAlleleInfo` (`models/classification_variant_info_models.py`)

Central tracker for a variant being imported:

| Field | Type | Description |
|-------|------|-------------|
| `imported_genome_build` | FK → GenomeBuild | Build the HGVS/coord was provided in |
| `imported_c_hgvs` | CharField | HGVS as imported |
| `imported_g_hgvs` | CharField | Genomic HGVS as imported |
| `status` | CharField | Matching status (Unresolved/Matched/Failed) |
| `grch37` | FK → ResolvedVariantInfo | Resolved GRCh37 info |
| `grch38` | FK → ResolvedVariantInfo | Resolved GRCh38 info |
| `allele` | FK → Allele | Matched allele (once resolved) |

### `ResolvedVariantInfo` (`models/classification_variant_info_models.py`)

Per-build resolution cache: HGVS → Variant mapping with validation status.

### `ClassificationImportRun` (`models/classification_import_run.py`)

Tracks import statistics per run: rows processed, created, updated, errors, warnings.

### `BulkClassificationInserter` (`classification_inserter.py` / `classification_patcher.py`)

- Accepts JSON payload conforming to EvidenceKey schema
- Creates or updates Classification records
- Validates field values against EvidenceKey definitions
- Handles publication at specified share level
- Returns validation responses per-record

---

## Export System

### CSV Export

Full or filtered export of classifications as CSV. Configurable columns.

### REDCap Export

Generates REDCap data dictionary + export format for research data capture.

### ClinVar Export

See ClinVar Integration section above.

### API Export (`ClassificationApiExportView`)

REST endpoint returning classifications in versioned JSON format (v1, v2, v3).

---

## API Endpoints

The classification REST API is versioned (v1, v2, v3) for backward compatibility.

### Key Endpoints

```
GET  /api/classifications/v1/<classification_id>   — Retrieve single classification
POST /api/classifications/v1/                      — Create classification
PUT  /api/classifications/v1/<classification_id>   — Update classification
GET  /api/classifications/v2/export/               — Bulk export with filtering
```

### Serialization

`serializers.py` provides DRF serializers for:
- Classification detail (full evidence)
- ClassificationModification (versioned snapshot)
- EvidenceKey definitions

---

## Celery Tasks

| Task | File | Description |
|------|------|-------------|
| `process_classification_import_task` | `classification_import_task.py` | Main import task for a `ClassificationImport` batch |
| `ClassificationImportMapInsertTask` | `classification_import_map_and_insert_task.py` | Runs omni_importer mapping then insert |
| `ClassificationImportProcessVariantsTask` | `classification_import_process_variants_task.py` | Processes VCF-inserted variants, links to classifications, triggers liftover |
| `ClassificationCandidateSearchTask` | `classification_candidate_search_tasks.py` | Finds classifications needing evidence updates or cross-sample review |

---

## Auto-Population System (`autopopulate_evidence_keys/`)

When a classification is created or a variant is resolved, evidence fields are automatically populated from:

- **Variant annotation** (VEP output): gnomAD AF, SpliceAI, SIFT, PolyPhen, consequence, etc.
- **ClinVar data**: ClinVar clinical significance, star rating, condition
- **Gene data**: Gene symbol, canonical transcript
- **ClinGen**: ClinGen allele ID

Configuration (`ClassificationEvidenceUpdateForm`) lets labs set thresholds:
- Minimum gnomAD AF to auto-apply BA1
- Minimum ClinVar star count to use
- Whether to overwrite existing values

Auto-population logic lives in `autopopulate_evidence_keys/` with one module per source.

---

## Grouping Models

### `ClassificationGrouping` (`models/classification_grouping.py`)

Database-backed grouping of classifications for the same variant/allele. Used for:
- Showing "all classifications for this allele" across labs
- Computing consensus significance
- Detecting discordance

### `AlleleGrouping` (`models/classification_grouping.py`)

Groups by allele origin (germline vs somatic). A single allele can have both germline and somatic classifications.

### `ClassificationLabSummary` (`models/classification_lab_summaries.py`)

Cached per-lab statistics:
- Total classifications by significance
- Classification trends over time
- Used for dashboard graphs

### `DiscordanceLabSummary` (`models/discordance_lab_summaries.py`)

Cached discordance statistics per lab for the dashboard.

---

## Admin Interfaces

All major models have Django admin classes for operational management:

- `ClassificationAdmin` — Create/edit/delete classifications
- `EvidenceKeyAdmin` — Define evidence schema
- `ClinicalContextAdmin` — Manage discordance contexts
- `DiscordanceReportAdmin` / `DiscordanceReportTriageAdmin` — Manage discordances
- `ClinVarExportAdmin` / `ClinVarAlleleAdmin` — Manage ClinVar submissions
- `ConditionTextAdmin` — Manage condition matching
- `ImportedAlleleInfoAdmin` / `ImportedAlleleInfoValidationAdmin` — Debug imports
- `ClassificationImportRunAdmin` — Track import runs
- `ClassificationReportTemplateAdmin` — Manage report templates

---

## Views Overview

### Core CRUD Views (`views/views.py`)

- Classification list with filtering
- Create classification (form or API)
- View individual classification (full evidence, history, diffs)
- History and version diff views
- Withdraw classification

### Dashboard (`views/classification_dashboard_view.py`)

- Lab statistics and trends
- Classification counts by significance
- Discordance overview
- VUS accumulation graphs

### Export (`views/classification_export_view.py`)

- CSV export with filter options
- REDCap data dictionary download
- API export endpoint

### ClinVar (`views/clinvar_export_view.py`)

- ClinVar export summary and detail
- Batch submission management
- Legacy ClinVar matching UI

### Discordance (`views/discordance_report_views.py`)

- List all discordance reports
- View individual report with classifications
- Triage/resolution UI for lab members
- Export discordance data

### Condition Matching (`views/condition_matching_view.py`)

- UI for matching condition text to ontology
- Testing matching algorithms
- Handling obsolete ontology terms

### Overlaps (`views/classification_overlaps_view.py`)

- Shows classifications shared between labs for same variants
- Lab agreement/disagreement analysis
- Patient overlap detection

### Candidate Search (`views/classification_candidate_search_view.py`)

- Finds classifications that may need review
- Searches for same variants classified differently in other samples

### Evidence Keys (`views/evidence_keys_view.py`)

- Browse all available evidence keys
- Lab-specific configuration

---

## URL Patterns (Key Groups)

```
/classification/activity/              — Activity logs
/classification/classifications/       — Classification list/search
/classification/create/                — Create new classification
/classification/<id>/                  — View classification detail
/classification/<id>/history/          — Version history
/classification/<id>/diff/             — Version diff
/classification/export/                — Data exports
/classification/import/                — Import tool
/classification/evidence_keys/         — Evidence key browser
/classification/clinvar_export/        — ClinVar export management
/classification/clinvar_match/         — Legacy ClinVar matching
/classification/condition_matching/    — Condition text matching
/classification/discordance_report/    — Discordance tracking
/classification/overlaps/             — Lab overlap analysis
/classification/dashboard/             — Analytics dashboard
/classification/groupings/            — Classification groupings
/classification/api/v1/               — REST API v1
/classification/api/v2/               — REST API v2
/classification/api/v3/               — REST API v3
```

---

## Forms

| Form | Purpose |
|------|---------|
| `EvidenceKeyForm` | Select evidence key + provide value |
| `ClassificationAlleleOriginForm` | Choose germline/somatic/other |
| `ClinicalSignificanceForm` | Filter by significance checkboxes |
| `ClassificationEvidenceUpdateForm` | Configure auto-population thresholds (gnomAD AF, ClinVar stars, etc.) |

---

## Templates (Key)

```
templates/classification/
├── classification_view.html          — Full classification detail view
├── classification_diff.html          — Version diff view
├── classification_history.html       — Version history list
├── classifications_legacy.html       — Classification list table
├── classification_export.html        — Export tool UI
├── classification_import_tool.html   — Import interface
├── evidence_keys.html                — Evidence key reference
├── classification_dashboard*.html    — Analytics dashboards
├── classification_groupings.html     — Grouped classifications
├── clinvar_export.html               — ClinVar export management
├── clinvar_export_detail.html        — ClinVar submission detail
├── clinvar_match_detail.html         — Legacy ClinVar matching
├── discordance_reports.html          — Discordance report list
├── discordance_report.html           — Report detail
├── discordance_report_triage_detail.html — Triage interface
├── condition_matching.html           — Condition text matching UI
├── imported_allele_info.html         — Import allele resolution detail
├── uploaded_classifications_unmapped.html — Upload status
└── emails/
    ├── classification_summary_email.html — HTML summary email
    └── classification_summary_email.txt  — Plain text summary email
```

---

## Key Enums

### `ClinicalSignificance` (`enums/classification_enums.py`)

| Code | Label |
|------|-------|
| `B` | Benign |
| `LB` | Likely Benign |
| `VUS` | Variant of Uncertain Significance |
| `LP` | Likely Pathogenic |
| `P` | Pathogenic |

### `SomaticClinicalSignificance`

AMP Tier system (Tier I–IV with levels A–D).

### `AlleleOriginBucket`

`Germline` / `Somatic` / `Unknown`

### `DiscordanceLevel`

```
no_entries → single_submission → concordant → discordant_minor → discordant_major
```

### `DiscordanceReportResolution`

`Ongoing` / `Concordant` / `Continued Discordance`

### `ShareLevel`

`user` → `lab` → `institution` → `logged_in_users` → `public`

---

## Signal Handlers (`signals/`)

Django signals trigger cascading updates:

- On `Classification` publish → recalculate `ClinicalContext` discordance
- On `ClassificationModification` save → update `ClassificationLabSummary` cache
- On allele resolved → link classification to allele, trigger liftover
- On condition text match updated → refresh related classifications

---

## Key Architectural Patterns

1. **Versioned Immutability** — `ClassificationModification` provides a complete, auditable history. The current `Classification` is a mutable pointer; history records are immutable.

2. **Evidence as Flexible JSON** — All classification data stored as key/value pairs in a JSON blob, governed by `EvidenceKey` definitions. This allows labs to have different fields without schema changes.

3. **Multi-level Sharing** — `ShareLevel` controls who sees what. Fields also have `max_share_level` for fine-grained control (e.g., patient data never shown to public even on public classifications).

4. **Async Variant Resolution** — Variant matching (HGVS → database variant) is async via Celery. `ImportedAlleleInfo` tracks pending/resolved status.

5. **Lab Collaboration via Clinical Context** — `ClinicalContext` aggregates all classifications for an allele; discordance is computed automatically on any change.

6. **Ontology Standardization** — Condition text entered by labs is mapped to standard ontology terms (MONDO, HPO, OMIM) for cross-lab comparison and ClinVar submission.

7. **Auto-population from Annotation** — Evidence fields are pre-filled from VEP annotation data, reducing manual entry for computational criteria (gnomAD AF, SpliceAI, etc.).

8. **REST API Versioning** — Three API versions (v1, v2, v3) maintained for backward compatibility with lab integrations and Shariant.
