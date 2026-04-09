# Patients App Research

## Purpose

The `patients` Django app manages patient demographics, links patients to samples and specimens, supports phenotype text matching against ontologies, and provides CSV-based import/export workflows. It also maintains a full audit trail of changes made to patient records.

---

## Key Models

### Patient

The central model representing an individual patient.

**Fields:**
- `family_code` - Family identifier string
- `first_name` - Patient first name
- `last_name` - Patient last name
- `date_of_birth` - Date of birth
- `date_of_death` - Date of death (nullable)
- `sex` - Enum: M (Male) / F (Female) / U (Unknown)
- `phenotype` - Free-text phenotype description (from `HasPhenotypeDescriptionMixin`)
- `affected` - Boolean indicating affected status
- `consanguineous` - Boolean indicating consanguinity
- `research_consent` - Boolean indicating consent for research use
- `medicare` - Medicare identifier string
- `billing_details` - Billing information
- Address fields - Street, city, state, postcode, country
- `telephone` - Contact phone number
- `fake_data` - FK to a fake data source (used for test/demo data generation)
- `_deceased` - Internal boolean backing the `deceased` property
- `external_pk` - OneToOneField to ExternalPK (nullable; links to external system ID)

**Mixins/Base classes:**
- `GuardianPermissionsMixin` - Object-level permissions via django-guardian
- `HasPhenotypeDescriptionMixin` - Adds `phenotype` text field and phenotype matching support
- `ExternallyManagedModel` - Supports records whose primary source of truth is an external system (LIMS etc.)
- `PreviewModelMixin` - Provides a `preview` property for lightweight display cards

**Methods:**
- `match()` - Finds a matching patient by name, date of birth, and sex. Returns a match type of EXACT (all fields match) or PARTIAL (some fields match). Used during CSV import to avoid creating duplicates.
- `get_samples()` - Returns all samples linked to this patient via the Specimen relationship chain.

**Properties:**
- `name` - Full name string (concatenated first + last)
- `deceased` - Boolean; derived from `_deceased` and `date_of_death`
- `age` - Computed age from `date_of_birth`
- `code` - Short display code for the patient
- `preview` - PreviewData instance for use in UI cards and autocomplete results

---

### Specimen

Represents a physical biological specimen collected from a patient.

**Primary key:** `reference_id` (string identifier, not auto-generated integer)

**Fields:**
- `description` - Free-text description of the specimen
- `collected_by` - Name or identifier of the collector
- `patient` - FK to Patient
- `tissue` - FK to Tissue lookup table
- `collection_date` - Date specimen was collected
- `received_date` - Date specimen was received by the lab
- `mutation_type` - Enum: G (Germline) / S (Somatic)
- `nucleic_acid_source` - Enum: D (DNA) / R (RNA)
- `_age_at_collection_date` - Stored age at time of collection (computed/cached)

---

### Tissue

Simple lookup/reference table for tissue types.

**Fields:**
- `name` - Tissue name
- `description` - Description of the tissue type

---

### PatientPopulation

Links a patient to one or more population groups (many-to-one with Patient).

Allows a single patient to be associated with multiple population ancestries or ethnic groups, which is relevant for population-stratified variant interpretation.

---

### PatientAttachment

File attachments associated with a patient record.

**Behavior:**
- Automatically generates image thumbnails for image attachments
- Uses `PrivateUploadStorage` to ensure files are not publicly accessible and are served through Django's access-controlled views

---

### Clinician

Represents a medical clinician who may be associated with patients or cases.

**Fields:**
- `email` - Contact email
- `first_name` - First name
- `last_name` - Last name
- `title` - Professional title (Dr, Prof, etc.)
- `specialty` - Clinical specialty
- `phone` - Phone number
- `user` - OneToOneField to Django User (nullable; not all clinicians have system accounts)

**Methods:**
- `match()` - Finds a clinician by a name string, supporting fuzzy/partial matching
- `cleaned_get_or_create()` - Parses a name string (handling titles, initials, etc.) and gets or creates the corresponding Clinician record
- `user_is_clinician()` - Class/static method checking whether a given User has an associated Clinician record

---

### ExternalPK

Stores an identifier assigned by an external system (such as a LIMS) for a patient or case.

**Fields:**
- `code` - The external identifier value
- `external_type` - The type/namespace of the external ID
- `external_manager` - FK to ExternalModelManager

**Unique constraint:** `(code, external_type, external_manager)`

---

### ExternalModelManager

A registry entry representing an external system that can own or manage records.

**Fields:**
- `name` - Display name of the external system
- `can_modify` - Boolean; whether VariantGrid may modify records owned by this system
- `modifications_sent_to_external_system` - Boolean; whether changes are pushed back to the external system

---

### PatientImport

The parent record for a CSV import event. Kept permanently for audit trail purposes even after the import is complete.

Each import produces one PatientImport, which in turn contains multiple PatientRecord rows (one per CSV row).

---

### PatientRecord

Represents a single row from a CSV import, in an intermediate state before being committed to the Patient/Specimen models.

**Fields:**
- `patient_records` - FK to PatientImport
- `record_id` - Row identifier within the import
- `valid` - Boolean indicating whether this row passed validation
- `validation_message` - Text describing any validation errors
- `sample` - FK to snpdb.Sample (if matched)
- `patient` - FK to Patient (if matched or created)
- `patient_match` - Enum: EXACT / PARTIAL (result of Patient.match())
- `specimen` - FK to Specimen (if matched or created)
- CSV data fields - Raw values from each CSV column, stored for reference

**Key method:** `process_patient_records()` - Iterates PatientRecord rows, calls `Patient.match()` to find or create patients, assigns samples, creates Specimens, and records a PatientModification for each change made.

---

### PatientModification

Audit trail entry recording a change made to a patient record.

**Fields:**
- `patient` - FK to Patient
- `user` - FK to Django User (who made the change)
- `date` - Timestamp of the modification
- `description` - Human-readable description of what changed
- `origin` - Enum: UPLOADED_CSV / MANUAL_VG_GUI / EXTERNAL_DATABASE / INTERNAL_VG
- `patient_import` - FK to PatientImport (nullable; set for CSV-originated changes)

---

## Patient to Sample Link

The chain connecting a Patient to sequencing data follows this structure:

```
Patient
  └── Specimen (patient FK; one patient has many specimens)
        └── snpdb.Sample (linked to Specimen)
```

During CSV import, `PatientRecord.process_patient_records()` orchestrates the full workflow:
1. Calls `Patient.match()` to find an existing patient or determine one should be created
2. Creates or updates the Patient record
3. Creates or updates the Specimen record (linked to the patient)
4. Links the Sample to the Specimen
5. Creates a PatientModification entry for the audit trail

---

## Phenotype Matching

### HasPhenotypeDescriptionMixin

Adding this mixin to a model provides a `phenotype` free-text field and hooks into the platform-wide phenotype matching pipeline.

### bulk_patient_phenotype_matching()

A management-level function that processes phenotype text for all patients (or a specified subset):
1. Parses free-text phenotype descriptions
2. Matches candidate HPO/MONDO terms against the ontology
3. Creates `TextPhenotypeMatch` entries associating patients with ontology terms
4. Flags matches for curator review via an approval workflow

### Approval Workflow

Curators review proposed TextPhenotypeMatch entries and approve or reject them. This ensures phenotype-ontology mappings are clinically validated before being used in analysis.

---

## URL Patterns

| URL | Purpose |
|-----|---------|
| `/patients/` | Patient grid/list view |
| `/patient/create/` | Create new patient |
| `/patient/<pk>/` | View patient detail |
| `/patient/<pk>/contact/` | Patient contact information |
| `/patient/<pk>/specimens/` | Patient specimens list |
| `/patient/<pk>/genes/` | Genes associated with patient variants |
| `/patient/<pk>/modifications/` | Patient modification audit trail |
| `/patient_imports/` | List of CSV import events |
| `/patient_import/<pk>/` | Detail view for a specific import |
| `/patient_import/<pk>/record/<record_pk>/` | Detail for a single CSV row record |
| `/patient/<pk>/file/upload/` | Upload a patient attachment |
| `/patient/<pk>/file/<file_pk>/delete/` | Delete a patient attachment |
| `/patient_term_matches/` | View phenotype term match results |
| `/bulk_patient_term/` | Trigger bulk phenotype matching |
| `/patient_term_approvals/` | Curator approval interface for term matches |

### Autocomplete Endpoints

| Endpoint | Purpose |
|----------|---------|
| `/autocomplete/Patient/` | Patient autocomplete |
| `/autocomplete/Specimen/` | Specimen autocomplete |
| `/autocomplete/Clinician/` | Clinician autocomplete |
