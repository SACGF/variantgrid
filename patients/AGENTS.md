# patients — agent notes
Owns: Patient and the material taken from them (Specimen, Extraction, SpecimenMeasure), externally managed records
(ExternalPK / ExternallyManagedModel), the patient records CSV import, extraction matching, and phenotype text matched
to HPO / OMIM / MONDO terms (PhenotypeDescription → TextPhenotypeSentence → TextPhenotype → TextPhenotypeMatch).
Start with:
- models/models_patient.py — Patient, Specimen, Extraction, imports and their audit rows
- models/models_phenotype.py — the phenotype text models, and patient_phenotype_terms for many patients in one query
- models/has_phenotype_description_mixin.py — what a Patient or Cohort gets for having phenotype text
- phenotype_matching.py — splitting text into sentences, matching and bulk matching; phenotype_matcher.py — the lookups
Patterns here:
- A PhenotypeDescription is owned through the one-to-one on the description, named after the owner's model
  (`patient`, `cohort`), so deleting the owner deletes it and its sentences. New text deletes the description and makes
  another, and the approval (`approved_by`) goes with it. An unowned description is a live preview
  (views_json.py:phenotypes_matches) or an SA Path request's, held by sapath's own link rows.
- TextPhenotype is one row per distinct sentence, shared by every description containing it and matched once.
  `match_version` null is awaiting matching; patients/phenotype_matching.py:requeue_sentences sets it back to null, and
  the old matches show until patients/phenotype_matching.py:bulk_patient_phenotype_matching replaces them.
- Ambiguous acronyms are one rule, patients/models/models_phenotype.py:TextPhenotypeMatch.is_ambiguous_acronym (and
  without_ambiguous_acronyms): matching doesn't save them, and every read path drops rows saved before the text
  joined the denylist. Read the denylist once per description or page and pass it down.
Gotchas:
- models_phenotype.py imports ontology, which imports genes and snpdb, which import patients.models - so
  patients/models/__init__.py star-imports only models_patient and patients/apps.py:PatientsConfig.import_models loads
  models_phenotype. Import the phenotype models from patients/models/models_phenotype.py, and keep models_patient.py
  and has_phenotype_description_mixin.py free of imports above snpdb.
- For the same reason the mixin can't import the matcher: process_phenotype_if_changed sends
  `phenotype_description_wanted_signal` and patients/signals/phenotype_description.py builds the description.
- patients/phenotype_matcher.py:get_ambiguous_acronym_denylist reads every ontology term and relation (~230MB) on a cache
  miss - 100s on a cold disk inside a page render. It is cached with no expiry and prebuilt on a new
  OntologyVersion (patients/tasks/ambiguous_acronym_denylist_task.py); a Redis flush means one slow rebuild.
- patients/phenotype_matcher.py:PhenotypeMatcher._get_special_case_lookups patches HPO gaps by name or ID, and goes
  stale as HPO adds and renames terms ("distal hypermobility" pointed at its opposite). After an HPO upgrade, compare
  each entry with the matcher's result without it; drop entries HPO now matches and repoint ones it contradicts.
- A change to phenotype lookups or matching logic bumps patients/phenotype_matcher.py `PHENOTYPE_MATCHER_VERSION` and
  adds a migration registering `match_patient_phenotypes --stale` as a ManualOperation, or deployments keep the old
  matches (sentences are matched once and cached). Each sentence points at the PhenotypeMatchVersion (matcher +
  OntologyVersion pair) it was matched with; stale and awaiting sentences show in `vg status`.
Tests:
- tests/test_phenotype_matching.py builds its ontology with ontology/tests/test_data_ontology.py:create_ontology_test_data
  and create_test_ontology_version. A PhenotypeMatcher takes seconds to build: make one in setUpTestData and pass it
  as `phenotype_matcher=` to `save`.
Deep reference: claude/research/patients.md · claude/maps/models.md#patients
