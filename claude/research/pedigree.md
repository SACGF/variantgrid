# Pedigree App Research

## Purpose

The `pedigree` Django app manages family pedigrees derived from PED files. It validates pedigree structures (family relationships, sex assignments, affected status), and links pedigree individuals to sequencing samples via cohort membership. The app also integrates with the Somalier tool for relatedness verification and provides pedigree chart rendering via the external Madeline2 tool.

---

## PED File Format

PED files are tab-delimited text files with six standard columns per row:

| Column | Field | Notes |
|--------|-------|-------|
| 1 | Family_ID | Groups individuals into families |
| 2 | Individual_ID | Unique identifier for this individual |
| 3 | Paternal_ID | Individual_ID of father; 0 if unknown |
| 4 | Maternal_ID | Individual_ID of mother; 0 if unknown |
| 5 | Sex | 1 = male, 2 = female, 0 = unknown |
| 6 | Affection_Status | 0 = unaffected, 1 = affected, -9 = unknown |

Each row represents one individual. Multiple families can appear in a single file, distinguished by their Family_ID values.

---

## Key Models

### PedFile

Represents an uploaded PED file.

**Fields:**
- `name` - Display name for this PED file
- `user` - FK to Django User (owner/uploader)
- `import_status` - Tracks the status of the file parsing process

**Mixins:**
- `GuardianPermissionsMixin` - Provides object-level permissions; only the owner and explicitly granted users can access the file

---

### PedFileFamily

Represents one family block within a PED file (all rows sharing the same Family_ID).

**Fields:**
- `ped_file` - FK to PedFile
- `name` - The Family_ID string from the PED file

**Methods:**
- `get_records_count()` - Returns the number of individuals in this family
- `validate()` - Runs structural validation and returns a list of error messages (empty list if valid)
- `is_valid()` - Returns True if `validate()` produces no errors
- `filter_for_user()` - Class method returning families accessible to a given user

**Properties:**
- `errors` - Cached list of validation error strings; populated by `validate()`

**Validation rules enforced:**
- At least 1 record must exist in the family
- Any individual listed as a father must have sex = male
- Any individual listed as a mother must have sex = female
- At least 1 individual must have affection status = affected

---

### PedFileRecord

Represents a single individual (one row) from a PED file.

**Fields:**
- `family` - FK to PedFileFamily
- `sample` - The Individual_ID string from the PED file
- `father` - Self-referential FK to PedFileRecord (nullable; the paternal individual)
- `mother` - Self-referential FK to PedFileRecord (nullable; the maternal individual)
- `sex` - Enum: M (male) / F (female) / U (unknown)
- `affection` - Nullable boolean (True = affected, False = unaffected, None = unknown)

Father and mother references are self-referential foreign keys resolved after all records for the family are created, since the referenced individuals must exist first.

---

### Pedigree

Connects a validated PED file family to a Cohort, linking pedigree individuals to actual sequencing samples stored in the system.

**Fields:**
- `user` - FK to Django User (owner)
- `name` - Display name for this pedigree
- `cohort` - FK to snpdb.Cohort (the set of samples in this pedigree)
- `ped_file_family` - FK to PedFileFamily (the family structure definition)

**Mixins/Base classes:**
- `GuardianPermissionsAutoInitialSaveMixin` - Automatically sets up object-level permissions on first save
- `PreviewModelMixin` - Provides a `preview` property for lightweight display cards and autocomplete

**Methods:**
- `get_samples(affected=None)` - Returns samples in this pedigree; if `affected` is provided, filters to only affected (True) or unaffected (False) samples
- `get_absolute_url()` - Returns the URL for the pedigree detail view

**Properties:**
- `genome_build` - Derived from the cohort's samples; the genome build (GRCh37/GRCh38) this pedigree's data is aligned to

---

### CohortSamplePedFileRecord

Junction table mapping each sample in the cohort to its position (individual) in the pedigree.

**Fields:**
- `pedigree` - FK to Pedigree
- `cohort_sample` - FK to snpdb.CohortSample (a specific sample's membership in a cohort)
- `ped_file_record` - FK to PedFileRecord (the individual in the pedigree this sample corresponds to)

This three-way mapping allows a single cohort to have its samples assigned to specific pedigree roles (proband, parent, sibling, etc.) as defined by the PED file structure.

---

### SomalierPedigreeRelate

Extends the `SomalierRelate` model to tie relatedness verification results to a specific pedigree.

**Relationship:** OneToOneField with Pedigree

**Method:**
- `write_ped_file()` - Writes a PED-formatted file for use as input to the Somalier tool, which verifies that the sample relatedness observed in the sequencing data matches the claimed family structure

---

## Parsing Flow

The PED file parsing process follows these steps:

1. Read the uploaded PED file line by line
2. Group rows by Family_ID to identify distinct family blocks
3. For each family block:
   a. Create a `PedFileFamily` record
   b. Create one `PedFileRecord` per row (individual)
   c. After all records for the family are created, resolve father/mother references by looking up the corresponding PedFileRecord by Individual_ID and setting the self-referential FK

The deferred resolution of parent references (step 3c) is necessary because a parent may be listed after a child in the file.

---

## Linking to Samples

Once a PED file is parsed and validated, a Pedigree can be constructed to connect the family structure to actual sequencing samples:

1. A `Cohort` is selected (or created) containing the relevant samples
2. A `Pedigree` is created linking the Cohort to a PedFileFamily
3. For each sample in the Cohort, a `CohortSamplePedFileRecord` is created associating that sample with the appropriate `PedFileRecord` (individual role)

### Auto-matching

`create_automatch_pedigree()` attempts to automatically match samples in a cohort to individuals in a PED file family by comparing sample names to Individual_ID values. Where names match, `CohortSamplePedFileRecord` entries are created automatically without manual assignment.

---

## Pedigree Chart Rendering

Pedigree charts are rendered using **Madeline2**, an external command-line tool that produces SVG pedigree diagrams from structured input data.

**Behavior:**
- The chart is generated by invoking Madeline2 as a subprocess with the pedigree data formatted appropriately
- Output SVGs are cached to avoid redundant re-generation on repeated requests
- Charts are served via the `pedigree_chart/<id>` URL

---

## URL Patterns

| URL | Purpose |
|-----|---------|
| `/pedigrees/` | List all pedigrees accessible to the current user |
| `/view_pedigree/<id>/` | Detail view for a specific pedigree |
| `/ped_files/` | List all uploaded PED files |
| `/view_ped_file/<id>/` | Detail view for a specific PED file, showing families and validation status |
| `/pedigree_chart/<id>/` | Serve the Madeline2-rendered SVG chart for a pedigree |
| `/create_pedigree_from_cohort_and_ped_file_family/<id>/` | Create a new Pedigree by associating a Cohort with a PedFileFamily |

### Autocomplete Endpoint

| Endpoint | Purpose |
|----------|---------|
| `/autocomplete/Pedigree/` | Pedigree autocomplete (for forms linking to a pedigree) |
