# #2131 Version phenotype matching so stale matches are rematched

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-09
Status: in progress

Phenotype matches are stored once per unique sentence (`annotation/models/models_phenotype_match.py:TextPhenotype`,
one row per sentence text, `processed` once matched) and never redone, so a matcher change or a new ontology import
only affects sentences first seen afterwards. The only rematch today is `match_patient_phenotypes --clear`, which deletes
every `TextPhenotype` and `PatientTextPhenotype` (losing `approved_by`) and cascades away the sentences under
`CohortTextPhenotype`, which `bulk_patient_phenotype_matching` never rebuilds because it only walks patients.

This plan records on each sentence what it was matched with, and rematches at the sentence level so patient and cohort
links and approvals survive. #2130 (fuzzy matching) lands on top of it and is the first version bump.

## Data

```python
# annotation/models/models_phenotype_match.py
class TextPhenotype(models.Model):
    text = models.TextField(primary_key=True)
    processed = models.BooleanField(default=False)
    # What the matches hanging off this sentence were produced with. Null = matched before #2131, so stale.
    matcher_version = models.IntegerField(null=True, blank=True)
    ontology_version = models.ForeignKey(OntologyVersion, null=True, blank=True, on_delete=SET_NULL)
```

```python
# annotation/phenotype_matcher.py
# Bump when a change to the lookups or matching logic would change the matches of an already matched sentence.
# TextPhenotype records the version a sentence was matched with; `match_patient_phenotypes --stale` redoes older ones.
# Every bump also registers that command as a ManualOperation in a new migration so deployments rematch (#2131).
# 1: special-case lookups repointed after the HPO review (2026-10)
PHENOTYPE_MATCHER_VERSION = 1
```

Migrations, both in `annotation/`: one adding the two fields, one one-off registering
`ManualOperation(task_id=ManualOperation.task_id_manage(["match_patient_phenotypes", "--stale"]), requires=["ontology-imported"], test=...)`
where the test is "any processed `TextPhenotype` exists" (at that point they all have a null `matcher_version`).
`ManualMigrationOutstanding.outstanding_task` treats one successful run as satisfying every registration before it,
so later bumps re-register the same task id in their own migration and a deployment that picks up several bumps at
once rematches once.

## Behaviour

1. **Stamping.** `annotation/phenotype_matcher.py:PhenotypeMatcher` captures `OntologyVersion.latest(validate=False)`
   at `__init__` (the ontology its lookups were built from). `annotation/phenotype_matching.py:_process_text_phenotype`
   sets `matcher_version = PHENOTYPE_MATCHER_VERSION` and `ontology_version` from the matcher when it sets `processed`.
   The sentences `create_phenotype_description` marks processed without matching (no alphanumerics) are stamped the
   same way, so they are never reported stale. A `TextPhenotype.mark_processed(ontology_version)` method is the one
   place that sets the three fields.
2. **Stale.** `TextPhenotype.stale_qs()` (static): `processed=True` and not
   (`matcher_version=PHENOTYPE_MATCHER_VERSION` and `ontology_version=OntologyVersion.latest(validate=False)`). With no
   ontology imported the current stamp is `None`, so the queryset compares against that.
3. **Requeue.** `annotation/phenotype_matching.py:requeue_sentences(text_phenotype_qs) -> int`, in one transaction:
   delete the `TextPhenotypeMatch` rows of those sentences, `update(processed=False)`, and invalidate the day-long
   `PhenotypeDescription.get_ontology_term_ids` memo for every description holding one of them (cache_memoize's
   `.invalidate(description)`; it applies `args_rewrite`). Patient and cohort links, `PhenotypeDescription`,
   `TextPhenotypeSentence` and `approved_by` are untouched. Phase 2 of `bulk_patient_phenotype_matching` already
   matches every `processed=False` sentence, cohort sentences included, so nothing else changes there.
4. **Command.** `patients/management/commands/match_patient_phenotypes.py`:
   - `--stale` requeues `TextPhenotype.stale_qs()`;
   - `--clear` requeues every processed sentence and stops deleting `TextPhenotype` / `PatientTextPhenotype` rows;
   - both then run `bulk_patient_phenotype_matching(cores=...)` as now, and the summary prints the requeued count
     alongside the before/after per-ontology counts.
5. **Reporting.** `library/vg/status.py` gains a `phenotype_sentences` section (`{"matched": M, "stale": N}`) rendered
   as one line naming `match_patient_phenotypes --stale` when N > 0. Deploys see the ManualOperation through
   `manual_outstanding` as usual.

## Docs

- `annotation/AGENTS.md` gotcha, next to the special-case one: a change to lookups or matching logic bumps
  `PHENOTYPE_MATCHER_VERSION` and adds a migration registering `match_patient_phenotypes --stale`; stale sentences
  show in `vg status`.
- `manual/__manual_readme.md` example and `claude/research/patients.md` "Phenotype text" paragraph now describe
  `--stale` (and that `--clear` rematches everything, keeping approvals).

## Tests (`annotation/tests/test_phenotype_matching.py`)

- A matched sentence carries the current `matcher_version` and the test `OntologyVersion`; an alphanumeric-free one
  is stamped too.
- `stale_qs`: a sentence stamped with an older matcher version, or another ontology version, is stale; a current one
  is not.
- Requeue then `bulk_patient_phenotype_matching`: the sentence's matches are rebuilt with the new stamp, and the
  patient's `PatientTextPhenotype` (with `approved_by`) and a cohort's `CohortTextPhenotype` still point at the same
  `PhenotypeDescription`.

## Scale

vg-test2 holds 58 sentences; SA Pathology holds tens of thousands. The requeue is three set-based statements, and the
rematch is the existing `--cores` pool, so the cost is the same as today's `--clear` without the re-registration walk.
