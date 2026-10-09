# patients — research notes

Verified against a96540a68 on 2026-09-25

The patients app holds the people a Sample was sequenced from and the material in between:
`patients/models/models_patient.py:Patient` → `patients/models/models_patient.py:Specimen` (one tissue at one timepoint) →
`patients/models/models_patient.py:Extraction` (the DNA or RNA arm off it) → `snpdb.Sample`. Around that sit the ways those rows
arrive (a CSV upload, a REST API for lab clients, references named by VCFs and seqauto records that are parked until
they resolve), phenotype text matched to ontology terms, and a per-patient audit trail. Fields, URLs, the task and the
command are in the generated maps ([models](../maps/models.md#patients), [urls](../maps/urls.md#patients),
[tasks](../maps/tasks.md#patients), [commands](../maps/commands.md), [signals](../maps/signals.md));
`claude/domain.md` has the vocabulary.

## Flows

### Identity: who a record is, and who may see it

Patient, Specimen and Extraction are `library/django_utils/guardian_permissions_mixin.py:GuardianPermissionsMixin`
models, but only Patient carries permission rows: Specimen, Extraction and `patients/models/models_patient.py:SpecimenMeasure`
point `get_permission_class` / `get_permission_object` at their patient, so granting a patient grants everything off it.
All three are `patients/models/models_patient.py:ExternallyManagedModel`s: an optional one-to-one `patients/models/models_patient.py:ExternalPK`
(unique on code, external_type, external_manager) plus a local reference named by `LOCAL_REFERENCE_FIELD`
(`patient_code` on Patient, `reference_id` on the other two). `ExternallyManagedModel.can_write` is false when the
`patients/models/models_patient.py:ExternalModelManager` says `can_modify=False`, and `Patient.filter_writable_for_user` applies the
same rule in SQL for grids. `ExternallyManagedModel.short_identifier` is the one identity rule for previews, search and
node chips.

A patient is shown by `patients/models/models_patient.py:Patient.display_identity`: the de-identified code alone when there is one
(showing a name beside it would re-identify the patient), else the name, else `Patient:<pk>`.
`Patient.display_identity_expression` is the same rule as a SQL `Case` so grids can sort, filter and export on it; the
two must change together.

### Patient → Sample, both ways round

A Sample reaches a patient either directly (`Sample.patient`, which the CSV import and the sample form set) or through
`Sample.extraction.specimen.patient`. `snpdb/models/models_vcf.py:Sample.save` fills an empty patient from the
extraction, but never replaces a different one - that disagreement is left for the sample form to report. So every
"samples of this patient" query is a union: `patients/models/models_patient.py:Patient.get_samples` on the model side and
`patients/sample_grouping.py:SOURCE_LEVELS` for the analysis grouping node. `patients/sample_grouping.py` is the single
implementation of which Samples a Sample / Extraction / Specimen / Patient reaches
(`get_sample_group`, `get_patient_sample_tree`), shared by the SampleNode above sample level, its canvas chips and the
node editor tree, so the node and the editor cannot disagree about how many VCFs are in play. Variant counts come from
cohort genotype stats rows, not the variant table (`get_sample_variant_counts`). The patient page runs
`snpdb/models/models_somalier.py:get_same_individual_relatedness` over `get_samples()` to warn of a likely sample swap.

### CSV patient records import

An uploaded CSV becomes `upload/tasks/import_patient_records_task.py:ImportPatientRecords` (web_workers), which makes a
`patients/models/models_patient.py:PatientImport` and a `patients/models/models_patient.py:PatientRecords` (one-to-one) and calls
`patients/import_records.py:import_patient_records`. The file must carry every `patients/models/models_patient.py:PatientColumns`
header. Each row runs `patients/import_records.py:process_record` in its own `transaction.atomic()`: a row that fails
rolls back to a single invalid `patients/models/models_patient.py:PatientRecord` (`create_failed_patient_record`) and the rest still
import; an unparseable value is a validation message on an otherwise imported row. There is no review-then-submit
step - rows are applied as they are read.

Per row: `match_sample` finds the Sample by id or name among samples the user can write;
`patients/models/models_patient.py:Patient.match` looks the patient up among those the user can see on last name (required), first
name, and DOB / sex where given - a blank DOB or Unknown sex on the stored patient still matches, and that is recorded
as PARTIAL rather than EXACT. Several matches fail the row, naming the patients (`#2037`). No
match creates the patient (`create_patient`, upper-cased names, permissions to the user's groups). Deceased / date of death, family code, patient code, affected and consanguineous are then updated;
phenotype text is appended under a dated "From Imported CSV" line rather than replaced. The specimen is found by
`(patient, reference_id)`; a reference_id already used by a *different* patient fails the row with
`specimen_patient_clash_message` (the lookup is unrestricted, so the other patient is only described if the user can
view it), and its fields are updated from the row by `patients/import_records.py:set_fields_from_row` - a blank column
is "no answer" and leaves the stored value alone. `patients/models/models_patient.py:Specimen.get_or_create_extraction` picks the extraction the row describes - the sample's
own extraction if it already has one, else one with the same nucleic acid, else an unnamed one - so a TSO 500 specimen's
DNA and RNA arms stay separate on re-import. `assign_patient_to_sample` / `assign_extraction_to_sample` write a
`patients/models/models_patient.py:PatientModification` for every change, hung off the PatientImport. Phenotype matching for touched
patients runs once in bulk at the end (`patients/phenotype_matching.py:bulk_patient_phenotype_matching`).

### REST API and external references

`patients/views_rest.py` exposes Patient / Specimen / Extraction viewsets (#1707) so a lab client can accession before
or alongside posting a run. `patients/serializers.py:ExternallyManagedModelSerializer` makes every create an upsert:
look the row up by its ExternalPK triple, else by its local reference (`_local_reference_q` - patient_code, or
`(parent, reference_id)`), within what the user can see; an ExternalPK that already belongs to a row of that model
the user cannot see is a 400 (`_external_pk_for`). Parents are named with
`patients/serializers.py:ExternalReferenceField`, a bare string (the local reference) or an object with
`reference_id` / `code` + `external_type` / `external_manager`, parsed by
`patients/external_references.py:ExternalReference.from_data`. An unknown external_manager is a 400 for non-superusers
(`PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY`), because creating one decides `can_modify` as a side effect of a typo.

`patients/external_references.py:resolve_reference` has three outcomes, not two: one match is MATCHED; more than one,
or the external and local halves each matching a different row, is NEEDS_ATTENTION; nothing is PENDING, because an
extraction legitimately arrives after the VCF that names it. The API treats anything but MATCHED as a 400
(`resolve_or_raise`) since a specimen has nowhere to live without its patient.

### Parked extraction claims and reconciliation

`patients/models/models_patient.py:ExtractionMatchMixin` (on `snpdb.Sample` and `seqauto.SequencingSample`) is a claim about which
Extraction a row belongs to that may not be resolvable yet: the reference as JSON, a `MatchStatus`, an error and the
date the claim was parked. `ExtractionMatchMixin.apply_extraction_match` never touches a row that already has its
extraction, and restarts the clock only when the reference changes. Claims are made by the VCF import
(`upload/vcf/vcf_import.py:assign_sample_extractions`, from upload metadata or, where
`PATIENT_EXTRACTION_SAMPLE_NAME_REGEX` is set, derived from the sample name with `derived=True`), by the seqauto link
serializer, and by the DRAGEN TSO500 imports (which claim a Specimen on `LibraryQC` / `DragenTSO500CombinedVariantOutput`
via seqauto's own specimen-claim mixin).

`patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions` (db_workers; hourly beat, and fired after
a Specimen or Extraction API create and after each DRAGEN import) re-resolves every parked Sample / SequencingSample claim as the
row's user (the VCF owner for a Sample), promotes PENDING older than `PATIENT_EXTRACTION_MATCH_PENDING_DAYS` to
NEEDS_ATTENTION, re-resolves the DRAGEN specimen claims, re-links LibraryQC and CombinedVariantOutput rows to runs,
sheets and samples that landed late, and finally copies a SequencingSample's extraction down to its Sample when the link
call arrived after the VCF. NEEDS_ATTENTION is what the health check
(`patients/signals/extraction_match_health_check.py:extraction_match_health_check`) and the unmatched extractions page
(`patients/views.py:unmatched_extractions`) count; PENDING is shown as context only.

### Phenotype text

Patient (and snpdb's Cohort) is a `patients/models/has_phenotype_description_mixin.py:HasPhenotypeDescriptionMixin`:
`Patient.save` pops the phenotype kwargs and matches the text to HPO / OMIM / MONDO terms unless
`check_patient_text_phenotype=False`. The text is a `patients/models/models_phenotype.py:PhenotypeDescription`
that points at its patient or cohort (a one-to-one each, at most one set), so deleting the owner - `Patient.merge`
included - deletes it; changed text replaces it, approval and all (#2135). It is split into sentences, each a
`patients/models/models_phenotype.py:TextPhenotype` shared by every description containing it and matched once.
The matcher lives in `patients/phenotype_matcher.py` and the splitting, matching and bulk path in
`patients/phenotype_matching.py`; the mixin reaches it through `phenotype_description_wanted_signal`
(`patients/signals/phenotype_description.py`), as the matcher needs ontology, which sits above snpdb.
`manage.py match_patient_phenotypes` registers every patient's sentences and matches those awaiting. Each matched
sentence points at the `patients/models/models_phenotype.py:PhenotypeMatchVersion` (the `PHENOTYPE_MATCHER_VERSION`
and OntologyVersion pair) it was matched with; `--stale` requeues those not on the current pair (`vg status` counts
them, and the awaiting ones) and `--clear` requeues every sentence. A requeued sentence keeps its matches until it is
rematched, and descriptions and approvals are untouched (#2131). Ambiguous acronyms are one rule,
`patients/models/models_phenotype.py:TextPhenotypeMatch.is_ambiguous_acronym`, for matching and every read path.
Curators approve a patient's matched text on the term approvals page (`patients/views.py:patient_term_approvals`,
`patients/views_json.py:approve_patient_term`, which records a PatientModification). The patients page graphs come
from `patients/templatetags/patient_graph_tags.py`.

## Why it is shaped this way

- **Specimen and Extraction are separate** (#1704) because one tumour block yields a DNA and an RNA arm, and one
  extraction can be sequenced more than once - so SequencingSample and Sample point at Extraction, not the reverse.
  `SpecimenMeasure` holds only what is measured on the material (pathologist tumour content, synced from SA Path's
  Mocha with `source_payload` kept); per-analysis numbers such as TMB / MSI belong to the seqauto
  CombinedVariantOutput because a re-analysis changes them (#1904).
- **Identifiers are either an ExternalPK or a local reference**, never forced into one canonical form, because some
  deployments have a LIMS to be managed by and some do not. A code only means something with its type, hence the
  triple.
- **Unresolved is parked, not rejected**, because feeds (sample sheet, VCF, accessioning, DRAGEN outputs) arrive in any
  order. The pending window separates the load race from a real mismatch that wants a human.
- **Permissions live on the patient only**, so there is one set of Guardian rows to grant and filter against.
- **Dates rather than ages**: age and deceased are derived from date of birth / death and collection date; the
  underscore fields (`_deceased`, `_age_at_collection_date`) exist only for when a date is not available, and
  `Patient.save` refuses both date_of_death and `_deceased` at once.

## History

- 2020: app bootstrapped around HPO / OMIM text matching of patient phenotypes; CSV patient records import with
  PatientModification auditing (2021 added specimen assignment auditing).
- 2025: CSV import handles rows with no specimen, and hides no-op modifications.
- Aug 2026: Specimen gets a surrogate pk (it used to be keyed on `reference_id`), Extraction is split out, tissue status
  moves onto the specimen (#1704, #1706); the Patient / Specimen / Extraction API and VCFs naming their extraction
  (#1707, #1559); SampleNode can group at patient / specimen / extraction level (`patients/sample_grouping.py`).
- Sep 2026: `patient_code` and "the code alone, never the name beside it" (#1860); a sample linked to an extraction gets
  its patient (#1875); SpecimenMeasure, then narrowed to pathology only (#1904); a failing CSV row rolls back to an
  invalid record instead of failing the upload; clearer specimen clash message (#422); somalier check on the patient
  page (#196).

## Traps

- **Two rules kept in step by hand**: `Patient.get_samples` and `SOURCE_LEVELS[PATIENT]` (both must union the direct
  and extraction paths), and `Patient.display_identity` with `display_identity_expression`.
- **`Sample.patient` can disagree with the extraction's patient** - `Sample.save` only fills an empty one. Query both
  paths, never just one.
- **Specimen reference_id is unique per patient in the database** (`unique_together` with patient) and the API, but the
  CSV import refuses a reference_id another patient already has. A bare reference can therefore match several
  specimens through the API / VCF path, which resolves to NEEDS_ATTENTION rather than picking one.
- **A CSV re-import updates the specimen from filled columns only** (`#2037`): a blank
  column leaves the stored value alone, so it cannot clear a description or reset `tissue_status`; a filled column
  replaces what is stored.
- **An API update does not fire reconciliation**; only a Specimen / Extraction create does, so a claim that an update
  would satisfy waits for the hourly beat.
- **`Clinician` is vestigial**: nothing references it beyond admin and an autocomplete.
