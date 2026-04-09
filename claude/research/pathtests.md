# Pathtests App Research

## Purpose

The `pathtests` Django app manages pathology (laboratory) testing workflows within VariantGrid. It covers test definitions with versioned gene lists, clinical cases, test orders, gene modification requests, and integration with external Laboratory Information Management Systems (LIMS). It acts as the bridge between clinical orders and genomic analysis, tracking a case from initial request through to reporting.

---

## Key Models

### PathologyTest

Defines a named pathology test (e.g., "Cardiomyopathy Panel"). The primary key is the test name itself (a string), making tests identifiable by human-readable name rather than an opaque integer ID.

**Fields:**
- `name` (PK) - Test name string; also the primary key
- `curator` - FK to Django User (nullable; the user responsible for maintaining this test)
- `deleted` - Boolean; soft-delete flag (deleted tests are hidden from normal views but preserved for historical orders)
- `empty_test` - Boolean; marks a placeholder test with no gene content yet

**Methods:**
- `can_write(user)` - Returns True if the given user has permission to modify this test
- `is_curator(user)` - Returns True if the given user is the assigned curator
- `get_active_test_version()` - Returns the currently active PathologyTestVersion for this test (via ActivePathologyTestVersion)
- `get_absolute_url()` - Returns the URL for this test's detail view

---

### PathologyTestVersion

A specific versioned snapshot of a PathologyTest, containing the gene list and enrichment kit assignment at a point in time. Once confirmed, a version is locked and cannot be changed; modifications require creating a new version.

**Fields:**
- `pathology_test` - FK to PathologyTest
- `version` - Integer version number (starts at 1, increments with each new version)
- `confirmed_date` - DateTime (nullable); when set, this version is locked and no further modifications are allowed
- `gene_list` - FK to a gene list model (defines which genes are tested)
- `enrichment_kit` - FK to an enrichment kit (the sequencing panel/capture kit used)

**Unique constraint:** `(pathology_test, version)`

**Properties:**
- `can_modify` - True if `confirmed_date` is None (version is still a draft)
- `can_confirm` - True if `can_modify` is True AND an enrichment kit has been assigned (both required before locking)
- `is_active_test` - True if this version is the one referenced by ActivePathologyTestVersion

**Methods:**
- `next_version()` - Creates a new PathologyTestVersion with version number incremented by 1, cloning the current gene list and enrichment kit as a starting point for further modification
- `set_as_active_test()` - Updates (or creates) the ActivePathologyTestVersion record to point to this version
- `replace_gene_list(new_gene_list)` - Replaces the gene list on this version (only valid if `can_modify` is True)

---

### ActivePathologyTestVersion

A one-to-one record per PathologyTest that points to the currently active PathologyTestVersion.

**Relationship:** OneToOneField on `pathology_test` → PathologyTestVersion

Acts as an indirection layer: the "active" version can be changed without altering any of the version records themselves.

---

### PathologyTestGeneModificationRequest

Tracks requests to add or remove a gene from a PathologyTestVersion. Provides a formal review workflow so gene changes are documented and approved before being applied.

**Fields:**
- `pathology_test_version` - FK to PathologyTestVersion (the version being modified)
- `outcome` - Enum: PENDING / ACCEPTED / REJECTED
- `operation` - Enum: ADD / REMOVE
- `gene_symbol` - FK to a gene symbol record
- `user` - FK to Django User (who submitted the request)
- `comments` - Free-text justification or notes for the request

---

### Case

Represents a clinical case — a single patient investigation that may involve one or more pathology test orders.

**Fields:**
- `name` - Case name or identifier
- `lead_scientist` - FK to Django User (nullable; the scientist leading the analysis)
- `result_required_date` - Date by which results are needed
- `report_date` - Date the final report was issued
- `patient` - FK to patients.Patient
- `details` - Free-text clinical details/notes
- `status` - Enum: OPEN / NO_TEST / CLOSED_SOLVED / CLOSED_UNSOLVED
- `workflow_status` - Enum tracking lab pipeline stage (see below)
- `investigation_type` - Enum: SINGLE_SAMPLE / TRIO / COHORT
- `external_pk` - OneToOneField to patients.ExternalPK (nullable; the LIMS case identifier)

**Mixins/Base classes:**
- `PreviewModelMixin` - Provides a `preview` property for display cards
- `ExternallyManagedModel` - Supports cases whose primary source of truth is an external LIMS

**Properties:**
- `is_open` - True if status is OPEN
- `is_closed` - True if status is CLOSED_SOLVED or CLOSED_UNSOLVED
- `code` - Short display code for the case

### Workflow Status Progression

```
NA
  └── SAMPLE_PROCESSING
        └── LIBRARY_PREP_COMPLETE
              └── SEQUENCING_COMPLETE
                    └── VCF_READY
                          └── ANALYSIS_COMPLETE
                                └── REPORTING
```

This progression mirrors the physical laboratory workflow from sample receipt through to final analysis reporting.

---

### CaseClinician

Associates one or more clinicians with a case.

**Fields:**
- `case` - FK to Case
- `clinician` - FK to patients.Clinician
- `specified_on_clinical_grounds` - Boolean; True if this clinician is listed for clinical/diagnostic grounds (as opposed to administrative or referral)

---

### PathologyTestOrder

Represents a specific order for a pathology test within the context of a case. Links the clinical request to the sequencing workflow and eventually to a completed analysis.

**Fields:**
- `case` - FK to Case (nullable; an order may exist before a case is created)
- `pathology_test_version` - FK to PathologyTestVersion (nullable; which test version was ordered)
- `custom_gene_list` - FK to a custom gene list (nullable; for ad-hoc orders not using a standard test)
- `user` - FK to Django User (who placed the order)
- `started_library` - DateTime when library preparation began
- `finished_library` - DateTime when library preparation completed
- `started_sequencing` - DateTime when sequencing began
- `finished_sequencing` - DateTime when sequencing completed
- `order_completed` - DateTime when the entire order was marked complete
- `experiment` - FK to seqauto experiment record (nullable)
- `sequencing_run` - FK to seqauto sequencing run record (nullable)
- `external_pk` - OneToOneField to patients.ExternalPK (nullable; the LIMS order identifier)

---

### PathologyTestOrderSample

Links a PathologyTestOrder to the sequencing sample produced for it.

**Relationship:** OneToOneField on `order` → snpdb.Sample

This one-to-one constraint ensures each order maps to exactly one sample record once sequencing data is available.

---

## Test Versioning Workflow

The complete lifecycle for creating and updating a pathology test:

1. **Create PathologyTest** - Define the test by name and assign a curator
2. **Create PathologyTestVersion** - The first version (version=1) is created automatically as a draft
3. **Assign gene list and enrichment kit** - Both required before the version can be confirmed
4. **Confirm the version** - Sets `confirmed_date`, locking the version against further edits; this version can now be set as active
5. **Set as active test** - Calls `set_as_active_test()`, updating ActivePathologyTestVersion
6. **Request gene modifications** - PathologyTestGeneModificationRequest records track proposed ADD/REMOVE operations with curator review
7. **Create next version** - When updates are needed, call `next_version()` to clone the current version with an incremented version number; the clone starts as a draft
8. **Repeat from step 3** for the new version

At any point, only the version referenced by ActivePathologyTestVersion is used for new orders.

---

## LIMS Integration

Both Case and PathologyTestOrder support integration with an external LIMS via their `external_pk` fields (OneToOneField to patients.ExternalPK).

### Lookup by External ID

Dedicated URL patterns allow the LIMS to push a request or look up a record using its own identifier without knowing VariantGrid's internal primary keys:

- `view_external_case/<ext_pk>/` - Look up a Case by LIMS case ID
- `view_external_pathology_test_order/<ext_pk>/` - Look up a PathologyTestOrder by LIMS order ID

### Workflow Status Synchronisation

The `workflow_status` field on Case is designed to be updated by the LIMS as the sample progresses through the lab pipeline. The seven statuses (NA through REPORTING) map to discrete milestones that the LIMS can set to reflect the current lab state, giving scientists visibility within VariantGrid without needing to query the LIMS directly.

---

## URL Patterns

### Web Views

| URL | Purpose |
|-----|---------|
| `/pathology_tests/` | List all pathology tests |
| `/manage_pathology_tests/` | Curator management interface for all tests |
| `/view_pathology_test/<name>/` | Detail view for a specific test (by name PK) |
| `/view_pathology_test_version/<pk>/` | Detail view for a specific test version |
| `/modify_pathology_test_version/<pk>/` | Edit interface for a draft test version |
| `/cases/` | List all cases |
| `/view_case/<pk>/` | Detail view for a specific case |
| `/view_external_case/<ext_pk>/` | Look up case by LIMS external ID |
| `/view_pathology_test_order/<pk>/` | Detail view for a specific test order |
| `/clinician_cases/<id>/` | Cases associated with a specific clinician |
| `/my_cases/` | Cases assigned to the currently logged-in scientist |
| `/scientist_cases/<id>/` | Cases assigned to a specific scientist |

### REST API Endpoints

| URL | Purpose |
|-----|---------|
| `/api/view_pathology_test_version/<pk>/` | Retrieve test version details as JSON |
| `/api/view_latest_pathology_test_version/<name>/` | Retrieve the active test version for a test by name as JSON |
