# #1870 — Patient CSV names the DNA and RNA arms of one specimen

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

[#1870](https://github.com/SACGF/variantgrid/issues/1870), split out of the TSO 500 plan (SACGF/variantgrid_sapath#431)
and left over from #1704. The patient CSV (`patients/import_records.py`) can put two extractions on one specimen but
cannot *name* them, so it cannot describe the arms the API (#1707) and the sample-name regex already talk about -
`2600000001C` (DNA) and `2600000001B` (RNA) off block `2600000001`.

## The problem

A CSV row is one sample. Its `SPECIMEN_*` columns (`patients/models.py:PatientColumns`) carry the specimen plus one
`Specimen Nucleic acid source (DNA/RNA)` value, and `patients/models.py:Specimen.get_or_create_extraction` turns that
into the row's `Extraction`: the sample's current extraction if it is on this specimen, else the specimen's extraction
with that nucleic acid, else an unnamed one it can fill in, else a new one. Nothing in the file sets
`Extraction.reference_id` or `Extraction.extraction_date`, so:

- The arms a lab knows as `2600000001C` / `2600000001B` come out of a CSV import as two unnamed extractions
  ("DNA", "RNA"). A VCF whose sample name ends in `2600000001C` parks a claim for an extraction with that
  `reference_id` (`upload/vcf/vcf_import.py:_derive_extraction_reference`, `PATIENT_EXTRACTION_SAMPLE_NAME_REGEX`)
  that a CSV-made extraction can never satisfy.
- Two DNA extractions off one specimen (a re-extraction, `test_extraction_date_distinguishes_a_re_extraction` in
  `patients/tests/test_specimen_extraction.py`) are indistinguishable to the CSV: a second DNA row with no sample
  silently lands on `.first()`.
- A row with no sample (the last two rows of `upload/test_data/patient_upload.csv`) can describe one arm only.
- The export (`patients/views_json.py:get_patient_upload_csv`, `SAMPLE_QUERYSET_PATH`) writes the nucleic acid and
  nothing else about the extraction, so a downloaded-then-reuploaded file loses the reference and date.

The API already expresses all of this: `patients/serializers.py:ExtractionSerializer` posts `specimen` (a reference),
`reference_id`, `nucleic_acid_source`, `extraction_date`, upserting on `(specimen, reference_id)`
(`_local_reference_q`), and `test_both_arms_post_against_one_specimen` in `patients/tests/test_patient_api.py` is the
two-arm case. The CSV gets the same three fields and the same matching rule.

### Drift from the issue and the brief

- **"Round-trips one Extraction per Specimen"** is no longer quite true. Since #1704/#1745 the sample discriminates:
  `upload/test_data/patient_upload.csv` already has `ONC042_DNA` and `ONC042_RNA` rows under `SPEC-042`, and
  `test_one_specimen_two_arms` in `upload/tests/test_import_patient_records.py` proves two extractions on one
  specimen. What the CSV cannot do is name them, date them, or hold two arms without a sample per arm. This plan is
  about naming, and the sample rule stays.
- **There is no "TSO 500 loader that creates the DNA and RNA arms".** Nothing in this repo or `variantgrid_sapath/`
  creates `Extraction` rows from TSO 500 output. The arms are created by the API
  (`patients/views_rest.py:ExtractionViewSet`) or by hand (`patients/forms.py:ExtractionForm`); the VCF import only
  *derives* a reference from the sample name and parks it - `test_a_derived_reference_that_resolves_to_nothing_creates_nothing`
  in `patients/tests/test_extraction_matching.py:SampleNameFallbackTest` says so in its name. The MoCHA sync
  (`variantgrid_sapath/sapath/tasks/import_mocha.py`) caches `MochaSampleExtraction` rows and updates patients and
  specimen measures, never `Extraction`. So the CSV becomes the *third producer* of named extractions, and what it
  must agree with is the regex's `extraction` group (`2600000001C`), so parked VCF samples resolve after a CSV import.
  Today the CSV import never triggers `patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions`;
  the API's `perform_create` does. Fixed below.
- **`patients/import_records.py`'s module docstring** describes a review-then-SUBMIT flow ("Display what's going to
  happen with the records ... Then, click SUBMIT") that does not exist: `process_record` parses, matches and writes in
  one pass, and `validation_message` is recorded on the `PatientRecord` without stopping the row. The new validation
  errors follow the existing convention (`patients/import_records.py:match_sample`: the message is recorded, the thing
  that failed is left `None`, the rest of the row is applied). The docstring is corrected in the same change.
- **All 21 columns are required today** (`patients/import_records.py:import_patient_records` raises on any missing
  header). Backwards compatibility for the new columns therefore needs an optional-column notion, which there is not
  yet.
- `claude/research/patients.md` has no `Verified against` header; its PatientRecord section is a lead, not a fact.
  The patients app has no CLAUDE.md of its own.

## Data

### `patients/models.py:PatientRecord` - three columns and one FK

The row's raw extraction values sit beside the raw specimen values it already keeps, and the matched extraction
beside the matched specimen, so the import page can show what each row named and what it landed on.

```python
class PatientRecord(TimeStampedModel):
    ...
    specimen = models.ForeignKey(Specimen, null=True, related_name='related_specimen', on_delete=SET_NULL)
    specimen_match = models.CharField(max_length=1, choices=PatientRecordMatchType.choices, null=True)
    extraction = models.ForeignKey(Extraction, null=True, related_name='related_extraction', on_delete=SET_NULL)
    extraction_match = models.CharField(max_length=1, choices=PatientRecordMatchType.choices, null=True)
    ...
    specimen_tissue_status = models.CharField(max_length=1, choices=TissueStatus.choices, null=True)
    # Renamed from specimen_nucleic_acid_source - the value has lived on Extraction since #1704
    extraction_nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, null=True)
    extraction_reference_id = models.TextField(null=True)
    extraction_date = models.TextField(null=True)   # as typed, like specimen_collection_date
    specimen_age_at_collection_date = models.IntegerField(null=True, blank=True)
```

Migration patients/migrations/0019_patientrecord_extraction.py: `RenameField` for the nucleic acid column, `AddField`
x3. No `Extraction` change - `reference_id`, `nucleic_acid_source`, `extraction_date` and `unique_together
("specimen", "reference_id")` are what the CSV writes to.

### `patients/models.py:PatientColumns` - the column set

```python
class PatientColumns:
    ...
    SPECIMEN_TISSUE_STATUS = 'Specimen Tissue status (Reference/Affected/Unknown)'
    SPECIMEN_AGE_AT_COLLECTION_DATE = 'Age at collection date (mutually exclusive to date of birth)'
    EXTRACTION_REFERENCE_ID = 'Extraction Reference id'
    EXTRACTION_NUCLEIC_ACID_SOURCE = 'Extraction Nucleic acid source (DNA/RNA)'
    EXTRACTION_DATE = 'Extraction date'

    COLUMN_DETAILS = [..., (EXTRACTION_REFERENCE_ID, "String", "Names one extraction of the specimen, eg the DNA "
                            "and RNA arms of a block - one row per extraction, repeating the specimen columns. "
                            "Unique within the specimen. Blank: the row's sample, else the nucleic acid and "
                            "extraction date, say which extraction it is"),
                      (EXTRACTION_NUCLEIC_ACID_SOURCE, "Nucleic Acid Source", ""),
                      (EXTRACTION_DATE, "Date", "Tells two extractions of the same acid apart")]
    COLUMNS = [c[0] for c in COLUMN_DETAILS]                 # header order; extraction columns last
    OPTIONAL_COLUMNS = {EXTRACTION_REFERENCE_ID, EXTRACTION_DATE}
    # Header a file may still carry from before #1870 -> the column it means
    LEGACY_COLUMNS = {'Specimen Nucleic acid source (DNA/RNA)': EXTRACTION_NUCLEIC_ACID_SOURCE}

    SAMPLE_QUERYSET_PATH = {...,
                            EXTRACTION_REFERENCE_ID: "extraction__reference_id",
                            EXTRACTION_NUCLEIC_ACID_SOURCE: "extraction__nucleic_acid_source",
                            EXTRACTION_DATE: "extraction__extraction_date"}
```

## CSV column design

**One row per extraction, in the existing flat file, `EXTRACTION_*` columns after the `SPECIMEN_*` ones.** A row
already repeats the patient and specimen columns for each sample (the two `ONC042` rows), the export is one row per
sample, file type detection keys on one header (`PATIENT_LAST_NAME` in
`upload/import_task_factories/import_task_factories.py`), and pandas reads a single header. A separate extraction
section would break all four for no gain. Shape-wise the row is `ExtractionSerializer`'s payload flattened: the
specimen reference plus `reference_id`, `nucleic_acid_source`, `extraction_date`.

The nucleic acid column is renamed from `Specimen Nucleic acid source (DNA/RNA)` to
`Extraction Nucleic acid source (DNA/RNA)` because it sits beside `Extraction Reference id` and `Extraction date`
and has described the extraction since #1704; the old header is accepted through `LEGACY_COLUMNS`. The constant
`SPECIMEN_NUCLEIC_ACID_SOURCE` becomes `EXTRACTION_NUCLEIC_ACID_SOURCE` (`grep` shows its users are
`patients/import_records.py`, the two test modules and `patients/models.py` itself).

A row with a `Specimen Reference id` and blank extraction columns means what it means today. A row with an
`Extraction Reference id` and no `Specimen Reference id` is a validation error (an extraction has nowhere to live;
`patients/serializers.py:SpecimenSerializer` gives the API's version of the same rule).

### Backwards compatibility

`import_patient_records`: rename `LEGACY_COLUMNS` headers on the DataFrame first, then require
`set(COLUMNS) - OPTIONAL_COLUMNS`, then add any absent optional column as `None` so `process_record` reads every
column the same way. A pre-#1870 file - the 21 old headers - imports exactly as today: one extraction per row, the
sample discriminating, `test_one_specimen_two_arms` still passing against a legacy copy (Tests, below).

## Matching rules on re-import

Specimen, unchanged: `(patient, reference_id)`, a reference held by another patient's specimen raises
(`TestProcessRecordSpecimenPatientMismatch` in `patients/tests/test_patients_edge_cases.py`).

Extraction - `Specimen.get_or_create_extraction(nucleic_acid_source=None, sample=None, reference_id=None,
extraction_date=None)`, in this order, first hit wins:

1. **`reference_id` given**: `extraction_set.filter(reference_id=...)` - the row's identity, the same key the API
   upserts on. Found: set `nucleic_acid_source` and `extraction_date` from the row where the row gives them (a blank
   cell leaves the stored value alone, the `set_fields_if_blank` convention of the specimen block; a differing value
   overwrites, as the API does). Not found: create with all three. A reference another specimen's extraction holds
   is not a clash - `reference_id` is unique per specimen (`test_same_reference_on_another_specimen_is_fine` in
   `patients/tests/test_specimen_extraction.py`).
2. **No `reference_id`, `sample` given and `sample.extraction` is on this specimen**: that extraction, acid updated
   in place - today's rule, kept so an existing file re-imports as before.
3. **No `reference_id`**: discriminate on what the row does give. Candidates are the specimen's extractions filtered
   by `nucleic_acid_source` when given and by `extraction_date__date` when given (dates are stored as midnight
   datetimes by `patients/import_records.py:parse_date`, the way specimen dates are). One candidate: it. Several:
   `AmbiguousExtraction` (a `ValueError` subclass in `patients/models.py`) naming them - "Specimen 2600000001 has 2
   DNA extractions (2600000001C, 2600000001D): name one with Extraction Reference id or Extraction date". None, and
   the row gave an acid: an extraction with no acid yet is the row's to fill in (today's rule). None, and the row gave
   nothing: whatever is there, if exactly one (today's rule). Otherwise create.

Step 1 beats step 2 deliberately: a file that names `2600000001B` on a sample currently linked to `2600000001C` is
correcting the link, and `patients/import_records.py:assign_extraction_to_sample` records the change as a
`PatientModification`. Step 3 matches named extractions too, so a legacy file (acid, sample, no reference) fills in
the sample link on arms the API accessioned - today's behaviour for `SPEC-042`. The one behaviour change in an
existing file is that two same-acid extractions with nothing to tell them apart now raise instead of picking
`.first()`.

`process_record` catches `AmbiguousExtraction` into `validation_messages`, leaves `extraction = None` (so no sample
link is made for that row), and carries on - the `match_sample` convention. `extraction_match` on the `PatientRecord`
is `CREATED` / `EXACT` like `specimen_match`.

`import_patient_records` calls `reconcile_pending_extractions.delay()` when any row created an extraction, so a VCF
imported before its CSV settles the way it does after an API post (`ExtractionViewSet.perform_create`).

## Validation errors

Recorded on the row; the rest of the row is applied. A new validation message, rather than a new `parse_*`, only
where it is new logic:

- `Extraction Reference id` without `Specimen Reference id` - "Extraction Reference id 'X' needs a Specimen
  Reference id"; no extraction, no sample extraction link.
- `AmbiguousExtraction` as above.
- Unparseable `Extraction date`: `parse_date`'s existing message. Nucleic acid outside DNA/RNA:
  `patients/import_records.py:parse_choice`'s existing message.
- Same reference twice in one file with different acids: last row wins, as two API posts would. Not an error.

## Export

`SAMPLE_QUERYSET_PATH` gains the two paths above; `get_patient_upload_csv` is otherwise unchanged (it writes
`PatientColumns.COLUMNS`, so the new columns and the renamed header come for free). An extraction with no sample is
not exported, as a specimen with no sample is not today - the export is a sample list. `example_upload_csv_all` ->
edit -> upload is the round trip: every row names its extraction, so step 1 matches and nothing is duplicated.

## Code changes

- `patients/models.py`: `PatientColumns` and `PatientRecord` as above; `Specimen.get_or_create_extraction` gains
  `reference_id` and `extraction_date` and the order in "Matching rules"; its docstring is the one place the rule is
  written down. `AmbiguousExtraction`.
- `patients/import_records.py`: header aliasing and optional columns in `import_patient_records`; `process_record`
  reads the three columns, passes them through, catches `AmbiguousExtraction`, writes the new `PatientRecord`
  fields; the reconcile call; the module docstring corrected to what the code does.
- `patients/templates/patients/view_patient_record.html`: three labelled rows under the specimen ones (rename the
  nucleic acid label). `patients/grids.py:PatientRecordColumns`: `extraction__reference_id` and `extraction_match`
  beside the specimen pair.
- `patients/views.py:import_patient_records_details` needs nothing - it renders `COLUMN_DETAILS`.
- `claude/research/patients.md` PatientRecord section: the three fields and the matching order, one paragraph.
- `scripts/vg map` after the model change.

## Tests

The matching order and the column handling are ours; pandas, `RenameField` and `TextChoices` are not.

`patients/tests/test_specimen_extraction.py`:

- `TestGetOrCreateExtraction` (the rule, no CSV): a reference finds its extraction and updates the acid; a
  reference not on the specimen creates one even when the specimen has an unnamed extraction of that acid; a
  reference beats the sample's current extraction; no reference with two same-acid extractions raises
  `AmbiguousExtraction`, and giving the date picks one; a re-extraction is created rather than matched when the
  date differs. The existing three cases stay.
- `TestPatientCSVExtractionRoundTrip` (the CSV path through `patients/import_records.py:process_record`): two
  sample-less rows `2600000001` / `2600000001C` DNA and `2600000001B` RNA make two named extractions, and importing
  them again makes none; a row naming `2600000001B` relinks the sample from `2600000001C` and a `PatientModification`
  says so; an extraction reference without a specimen reference is a validation message and no extraction; an
  ambiguous row is a validation message, no sample link, and the specimen is otherwise updated; a row under the
  legacy nucleic acid header (built with `_make_row` and the old key) behaves as the new one;
  `test_export_columns_round_trip` reads the two new paths back.

`upload/tests/test_import_patient_records.py` (the whole pipeline, `TestPatientUploadImport`): regenerate
`upload/test_data/patient_upload.csv` as #1745 did - the renamed header, the two new columns, and the `ONC-042`
rows become block `2600000001` with arms `2600000001C` (DNA, dated) and `2600000001B` (RNA, dated), plus a third
sample-less row for a re-extraction `2600000001D` (DNA, later date). `test_one_specimen_two_arms` asserts the
references and dates; a new `test_legacy_header_imports` writes the same file back to the old shape in
`setUpTestData` (pandas: drop the two columns, rename the header, `tempfile`) and runs it through the pipeline,
asserting the old outcome (two arms under the sample rule, the re-extraction row now a validation message since it
is a second unnamed DNA). That single test is the whole backwards-compatibility promise; the file it needs is
derived, so there is no second CSV to keep in step.

## Manual verification

On this box, as a user with a VCF whose sample names end in `..._2600000001C` / `..._2600000001B` (the
`upload/test_data/tso500` example, `PATIENT_EXTRACTION_SAMPLE_NAME_REGEX` set as the setting's comment shows):

1. `/patients/help/import_patient_records_details` lists the three new columns with their descriptions.
2. `/patients/example_upload_csv/empty` downloads the new header. Fill two rows: one patient, specimen
   `2600000001`, extractions `2600000001C` DNA and `2600000001B` RNA, no samples. Upload; the import page shows
   both rows valid with Extraction Match CREATED; the patient's Extractions tab shows both arms named and dated.
3. Import the TSO 500 DNA VCF (it will have parked a claim if imported before step 2, or resolves at import): the
   sample's extraction is `2600000001C`. If it was parked, the reconcile fired by the CSV import settled it -
   `manage.py vg inspect sample <id>`.
4. Re-upload the same CSV: no new extractions, Extraction Match MATCHED EXACT.
5. `/patients/example_upload_csv/all`, then upload the file unchanged: same count of extractions; then change one
   `Extraction date` and re-upload: the date changes in place.
6. Upload the pre-#1870 `upload/test_data/patient_upload.csv` from `git show 36e3f4e5b:upload/test_data/patient_upload.csv`:
   imports as it did, seven rows, two arms under `SPEC-042`.

## Decisions made

- **Flat rows, one per extraction, `EXTRACTION_*` columns**, not a second section - the API's payload flattened
  onto the row shape the file, the export and the file-type detection already have.
- **Rename the nucleic acid header** and accept the old one as an alias. Leaving `Specimen Nucleic acid source`
  beside `Extraction Reference id` would have the help page contradict itself; the alias is a one-line
  `df.rename` and a test.
- **The two new columns are optional; everything else stays required.** Making all columns optional was
  tempting and out of scope - the strict header check is what stops a half-edited spreadsheet from silently
  importing.
- **Reference beats sample beats discriminators.** The reference is the identity the API upserts on; the sample
  rule is what existing files rely on; the discriminators are for rows that have neither.
- **Ambiguity is a validation message, not a guess.** Today's `.first()` is silent data loss for a re-extraction;
  a message on the row is the existing way this importer says "tell me more".
- **Unnamed extraction date on the `PatientRecord` is text**, like the specimen dates there - it records what was
  typed, and `Extraction.extraction_date` holds the parsed value.
- **The CSV import fires `reconcile_pending_extractions`** because the point of naming arms is that parked VCF
  samples can find them.

## Definition of done

- `scripts/vg tests --explain` names `patients.tests.test_specimen_extraction`,
  `patients.tests.test_patients_edge_cases` and `upload.tests.test_import_patient_records`; all pass with
  `--keepdb`; the kept tests are those above.
- Migration 0019 applied on this box (ask first); `scripts/vg map` refreshed.
- `Specimen.get_or_create_extraction`'s docstring carries the matching order; `patients/import_records.py`'s
  module docstring describes the one-pass import; `claude/research/patients.md` updated;
  `scripts/vg docs check` passes.
- This plan's `Status:` updated.
