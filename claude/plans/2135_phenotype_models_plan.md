# #2135 Phenotype models: owned descriptions, int sentence PK, one home in `patients`

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-09; revised to match the implementation by Claude Opus 5.5 (claude-opus-5-5), 2026-10-09
Status: in progress

The phenotype text models (then annotation/models/models_phenotype_match.py) are checked against a local copy of a
production database (17,973 sentences, 38,840 matches, 89,352 descriptions, 75,244 of them owned by nothing).
Findings on top of the issue:

- The orphans come from deleting a patient (`Patient.merge` included) and from the pre-#2131 `--clear`, which deleted
  every `PatientTextPhenotype`. Editing phenotype text does **not** orphan: `process_phenotype_if_changed` deletes the
  description and the link row cascades with it.
- Putting a nullable `phenotype_description` FK on `Patient`, as the issue suggests, does not fix the orphans: Django
  cascades along reverse relations only, so deleting the patient would still leave the description. The FK goes the
  other way - the description points at its owner - and the cascade is Django's.
- `TextPhenotype.processed` is the same fact as `match_version is not null` once the #2131 rematch has run. One field.
- `PhenotypeDescription.status` is `''` on every row; `DescriptionProcessingStatus` has a table and no columns.
- Reading the matches: `OntologyVersion.latest()` is an uncached query, run once per match in
  `TextPhenotypeMatch.to_dict`; `TextPhenotypeSentence.get_results` adds the sentence offset to the match instances.
- `snpdb` imports `patients.models` at module load, and the phenotype models import `ontology`, which imports `genes`
  and `snpdb`. So `models_patient.py` and the mixin stay below `snpdb`, and the phenotype models sit above it: the
  owners are string FKs and the mixin reaches the matcher through a signal, with every import at the top of its file.
- SA Path's `sapath` app (sibling repo) holds descriptions through its own link rows
  (`SAPathologyRequestConditionPhenotype`, `SAPathologyRequestClinicalNotesPhenotype`), so "owned" also counts any
  other relation to the description.

This plan does all seven items of the issue in one PR. #2131 is merged; #2130 (PR #2133) was merged into the #2131
branch after that had merged, so its two commits come in with this branch.

## Data

The models package `patients/models/`: `models_patient.py` (the old `patients/models.py`), `models_phenotype.py`
(below), `has_phenotype_description_mixin.py` (moved from annotation), and an `__init__.py` that star-imports
`models_patient` only. `patients/apps.py:PatientsConfig.import_models` loads `models_phenotype`, and callers import
the phenotype models from it.

```python
# patients/models/models_phenotype.py

class PhenotypeMatchVersion(TimeStampedModel):        # unchanged
    matcher_version = models.IntegerField()
    ontology_version = models.ForeignKey(OntologyVersion, null=True, blank=True, on_delete=CASCADE)
    # UniqueConstraint(matcher_version, ontology_version, nulls_distinct=False)


class TextPhenotype(models.Model):
    """ One row per distinct sentence - the unit of matching, shared by every description containing it """
    # id = AutoField (was: text as primary key)
    text = models.TextField(unique=True)
    # Null = awaiting matching (replaces `processed`); stamped when matched, set back to null to requeue
    match_version = models.ForeignKey(PhenotypeMatchVersion, null=True, blank=True, on_delete=SET_NULL)


class PhenotypeDescription(models.Model):
    """ A patient's or cohort's phenotype text, split into sentences. Owned by its patient or its cohort, so it
        goes when they do; owned by neither only for the moment a live preview needs one """
    original_text = models.TextField()
    patient = models.OneToOneField("patients.Patient", null=True, blank=True, related_name="phenotype_description",
                                   on_delete=CASCADE)
    cohort = models.OneToOneField("snpdb.Cohort", null=True, blank=True, related_name="phenotype_description",
                                  on_delete=CASCADE)
    approved_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    # CheckConstraint "phenotype_description_one_owner": patient is null or cohort is null


class TextPhenotypeSentence(models.Model):             # unchanged
    phenotype_description = models.ForeignKey(PhenotypeDescription, on_delete=CASCADE)
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    sentence_offset = models.IntegerField()


class TextPhenotypeMatch(models.Model):                # unchanged
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    offset_start = models.IntegerField()
    offset_end = models.IntegerField()
```

Gone: `DescriptionProcessingStatus`, `PhenotypeDescription.status`, `TextPhenotype.processed`, `PatientTextPhenotype`,
`CohortTextPhenotype`. `Patient.phenotype` and `Cohort.phenotype` are unchanged; approval lives on the description, so
new text (a new description) needs approving again, exactly as the link row did.

`PatientPhenotypeTerms` (dataclass) keeps its members and stays in `models_phenotype.py`.

## Migrations

Each is its own migration, applied to the dev database on this box before the PR, with `makemigrations --check`
clean at the end and the `--keepdb` test database dropped so it is rebuilt.

1. `patients/migrations/0020_phenotype_models_from_annotation.py`: `SeparateDatabaseAndState`. Database: `RunSQL`
   renaming each of the eight tables (`annotation_<model>` → `patients_<model>`). State: `CreateModel` for the eight
   models exactly as annotation had them, depending on annotation's 0192 and the current `snpdb` and `ontology`
   migrations. Index, constraint and sequence names keep their `annotation_` prefix; the stale `django_content_type`
   rows are left for `remove_stale_contenttypes`.
2. `annotation/migrations/0194_phenotype_models_to_patients.py`: state-only `DeleteModel` for each, children first.
   sapath's `0022_phenotype_description_moved_to_patients` (state-only `AlterField` of its two link FKs) runs between
   the two (`run_before` this one).
3. `patients/migrations/0021_phenotype_description_owner_fields.py`: `AddField` `patient`, `cohort`, `approved_by`.
   Depends on annotation's 0194, so sapath's links are in the state.
4. `patients/migrations/0022_phenotype_description_owners_copied.py`: data only, ORM. A `Subquery` update copies
   owner and `approved_by` from the two link tables; then descriptions nothing holds - no owner and no row of any
   other relation in the migration state, sapath's links included - are deleted with their sentences (75,244 locally).
   Data steps get their own migration so their deferred FK checks run at its commit, not inside the next one's
   `ALTER TABLE`s.
5. `patients/migrations/0023_phenotype_description_owner_links_removed.py`: `DeleteModel` the two link tables,
   `AddConstraint` one-owner.
6. `patients/migrations/0024_textphenotype_rebuilt_with_integer_key.py`: `RemoveField` `PhenotypeDescription.status`,
   `DeleteModel` `DescriptionProcessingStatus`, then `DeleteModel` / `CreateModel` `TextPhenotype`,
   `TextPhenotypeSentence` and `TextPhenotypeMatch` with integer keys - nothing is copied, as every sentence awaited a
   rematch anyway (processed before #2131's stamping). A `match_patient_phenotypes --rebuild` ManualOperation splits
   every description again (`patients/phenotype_matching.py:register_unsplit_descriptions`) and matches the
   sentences; until it runs, patients show no terms, so it runs straight after `migrate`.
7. `CACHE_VERSION` in `variantgrid/settings/components/default_settings.py` goes up: the modules of these classes
   change, so anything Redis pickled under the old paths fails to unpickle.

## Behaviour

### Mixin (`patients/models/has_phenotype_description_mixin.py`)

The subclass contract is a `phenotype` TextField and being the target of one owner field on `PhenotypeDescription`;
the owner field name is `self._meta.model_name` (`"patient"`, `"cohort"`), so the two abstract hooks and the inline
imports go (sapath's `SAPathologyRequest` overrides `get_phenotype_description` and `process_phenotype_if_changed`
for its link rows). `get_phenotype_description()` returns the description or `None` (the reverse one-to-one accessor,
`RelatedObjectDoesNotExist` → `None`, so a second call is free). `get_ontology_term_ids()`, `get_gene_symbols()`,
`pop_kwargs`, `save_phenotype` and the `save(**kwargs)` names (`check_patient_text_phenotype`,
`phenotype_approval_user`, `phenotype_matcher`) keep their signatures. `process_phenotype_if_changed` deletes the
description when the text differs (sentences cascade, approval goes with it) and otherwise sends
`phenotype_description_wanted_signal`; `patients/signals/phenotype_description.py` receives it and calls
`create_phenotype_description(text, owner=…, approved_by=…, phenotype_matcher=…, defer_processing=…)`.

Import layering, all at the top of each file: `snpdb` → `models_patient` → mixin (nothing above `snpdb`);
`models_phenotype` → `models_patient`, `ontology`, `patients/phenotype_matcher.py`; `patients/phenotype_matching.py`
→ `models_phenotype`; the signal receiver → `phenotype_matching`. `snpdb/models/models_cohort.py`,
`models_vcf.py` and `models_somalier.py` import from `models_patient` and the mixin module directly.
`scripts/vg imports cycles` and `lint-imports` stay clean.

### Matching (`patients/phenotype_matching.py`, `patients/phenotype_matcher.py`, `patients/phenotype_tokenizer.py`)

Moved from `annotation/` with their tests. `create_phenotype_description` gains `owner` and `approved_by` and
creates the description with them; the live preview (`patients/views_json.py:phenotypes_matches`) still creates an
unowned one and deletes it. "Awaiting" is `match_version__isnull=True`: `_process_text_phenotype` deletes the
sentence's existing matches before `bulk_create` (so a rematch is idempotent and a half-finished run leaves no
duplicates), then stamps `match_version`; `requeue_sentences` is `update(match_version=None)` plus the memo
invalidation - old matches stay visible until the sentence is rematched. `TextPhenotype.stale_qs()` is "stamped and
not current"; `TextPhenotype.awaiting_qs()` is "unstamped"; `library/vg/status.py` reports `matched`, `stale` and
`awaiting`.

`bulk_patient_phenotype_matching(patients, cores=1)` takes the queryset; the default (`phenotype` set and free of
`settings.PATIENT_PHENOTYPE_EXCLUDE_STRING`) becomes `Patient.with_phenotype_text()` and is built by the command,
`patients/management/commands/match_patient_phenotypes.py` (`--clear` requeues `match_version__isnull=False`).

The denylist celery task and its signal move to `patients/tasks/` and `patients/signals/` (registered in
`PatientsConfig.ready`), queue `db_workers` unchanged.

### Ambiguous acronyms, one rule

`TextPhenotypeMatch.is_ambiguous_acronym(denylist)` - `match_text` lowercased with commas stripped, the normalisation
`filter_ambiguous_acronym_matches` uses today - and `without_ambiguous_acronyms(matches, denylist)` over any
iterable of matches, saved or not. Used by `_process_text_phenotype` (before saving), `get_ontology_term_ids`,
`patient_phenotype_terms` and `TextPhenotypeSentence.get_results` (a row matched before its text joined the denylist
shows as a warning, as the term paths already treat it). `get_ambiguous_acronym_denylist()` is read once per call at
the description or bulk level and passed down.

### Results

`PhenotypeDescription.get_results()` resolves `OntologyVersion.latest()`, the denylist and the compiled acronym
pattern once and hands them to each sentence; `TextPhenotypeSentence.get_results(...)` builds new dicts with the
sentence offset added, leaving the match instances alone; `TextPhenotypeMatch.to_dict(ontology_version)` takes the
version. The JSON keys the templates and `patient_phenotype.js` read (`accession`, `gene_symbols`, `match`,
`ontology_service`, `name`, `offset_start`, `offset_end`, `pk`, `ambiguous`, `ambiguous_alias`,
`ambiguous_alias_candidates`) are unchanged. `PhenotypeDescription.__str__` is the owner and the first 50 characters.

### Patient-side helpers

`patients_qs_for_ontology_term` becomes `Patient.for_ontology_term(user, term)` in `models_patient.py`, with
`PATIENT_ONTOLOGY_TERM_PATH`; `patient_phenotype_terms(patients)` and `patient_phenotypes_for_samples(user, samples)`
stay in `models_phenotype.py`, which imports `Patient`. The path constants shrink by one hop:
`PATIENT_ONTOLOGY_TERM_PATH = "phenotype_description__textphenotypesentence__text_phenotype__textphenotypematch__ontology_term"`,
`TPM_PATIENT_PATH = "text_phenotype__textphenotypesentence__phenotype_description__patient"`.

### Callers

Every importer of `annotation.models.models_phenotype_match`, `annotation.models.has_phenotype_description_mixin`,
`annotation.phenotype_matching` and `annotation.phenotype_matcher` moves to the `patients` paths (grep; the
`annotation/models/__init__.py` star import goes). Relation callers: `patients/views.py:patient_term_approvals`
filters `phenotype_description__isnull=False` / `phenotype_description__approved_by__isnull=True` and
`select_related("phenotype_description")`; `patients/views_json.py:approve_patient_term` sets
`patient.phenotype_description.approved_by`; `patients/templates/patients/patient_term_approvals.html`,
`analysis/views/views_wizard.py`, `mme/serializers/patient_profile.py` (and the `SimpleNamespace` in
`mme/tests/test_ontology_snake.py`) read `phenotype_description` off the patient. `patients/models/models_patient.py:Patient.merge`'s
comment names the description. `snpdb/tests/test_patient_phenotypes_page.py` looks for the
`patients_textphenotypematch` table. The comment at `variantgrid/settings/components/default_settings.py` line 457
names the new models.

## Docs

- New patients/AGENTS.md: the three phenotype gotchas from `annotation/AGENTS.md` (denylist cost, special-case
  lookups, matcher version bump), plus: a description is owned through the FK on the description, so deleting the
  owner cascades; an unowned description exists only during a preview.
- `claude/research/patients.md` "Phenotype text", `claude/domain.md` (Cohort paragraph), `analysis/AGENTS.md`,
  `manual/__manual_readme.md` point at the new paths. `scripts/vg docs check` and `scripts/vg map` run clean.

## Tests

The phenotype matching tests move to `patients/tests/test_phenotype_matching.py` with their patch targets; `processed`
assertions become `match_version` ones; `test_requeue_keeps_links_and_approvals` reads
`patient.phenotype_description.approved_by`. `snpdb/tests/test_cohort_phenotype.py` reads
`cohort.phenotype_description`. New, each covering a rule above:

- deleting a patient deletes its description and sentences; the shared `TextPhenotype` stays;
- a description with both owners is refused by the constraint;
- matching a sentence that already has matches leaves exactly one set;
- `is_ambiguous_acronym` on a saved match agrees with the write-time filter for a comma-bearing match text.

Migrations are checked on the dev database: 14,108 descriptions keep their patient and `approved_by`; zero are
unowned afterwards; sentence, match and `TextPhenotype` counts are unchanged; `makemigrations --check` is clean.
`vg page` on a patient page and the cohort page, `--queries` before and after, shows the match-reading query count
falling.

## Scale

Every migration is set-based SQL over at most ~100k rows (seconds). Matching cost is unchanged: sentences are still
deduplicated and matched once per `PhenotypeMatchVersion`.
