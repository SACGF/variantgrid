# Consistent, then mandatory, mode of inheritance at SA Path

Written by Claude Opus 5 (claude-opus-5), 2026-09-16
Status: draft

The mode-of-inheritance half of [sapath#363](https://github.com/SACGF/variantgrid_sapath/issues/363)
(Classifications - condition + MOI). The condition matching half went to
[sapath#449](https://github.com/SACGF/variantgrid_sapath/issues/449).

The 2025-10-23 triage agreed three steps: MOI is a defined list, the legacy SA Path values are converted to
it, and then Karin (SA Path) is asked to make it mandatory. This plan is those three steps. Every choice
below that someone might want to make differently is listed under Decisions.

---

## What is already there

- **The list is already defined.** `mode_of_inheritance` is a multi-select (`value_type='M'`) with 18
  snake_case options and `allow_custom_values=False`, seeded in `classification/migrations/0002_blat_keys.py`.
  No later migration changes its options (`classification/migrations/0030_clinvar_mappings.py` only adds
  ClinVar mappings, `classification/migrations/0177_evidence_key_copy_scope_and_allele_origin_values.py` only
  copy scope), and nothing in the sapath repo configures it. So there is no evidence key change.
- **Every write today normalises.** Form, API, file import (`classification/models/classification_inserter.py:BulkClassificationInserter`)
  and the Shariant download (`sync/shariant/variant_grid_download.py`) all go through
  `classification/models/classification.py:Classification.patch_value` -> `Classification.process_entry`, which
  splits a string on `,` and turns each part into its option key with `Classification.process_option_values`;
  a part that matches nothing is an `invalid_value` error, which blocks submitting. So no new dirty value can
  arrive. The dirty rows are legacy data written before that path existed (VG3).
- **Matching** (`classification/models/classification.py:Classification.match_option`) compares
  case-insensitively against the key and `classification/models/evidence_key.py:EvidenceKey.pretty_label_from_string`
  (`x_linked_recessive` -> "X-linked recessive"), and again with `-` and spaces turned into `_`.

### The dirty values from the issue, checked

Run through `Classification.match_option` against the live key on vg-test2 (2026-09-16). **All 29 map**, so
no aliases are needed:

| Stored at SA Path | Count | Maps to |
|---|---|---|
| `Autosomal recessive` | 1039 | `autosomal_recessive` |
| `Autosomal dominant` | 495 | `autosomal_dominant` |
| `X-linked`, `x-linked` | 27, 1 | `x_linked` |
| `X-linked recessive` | 19 | `x_linked_recessive` |
| `X-linked dominant` | 10 | `x_linked_dominant` |
| `Isolated cases` | 6 | `isolated_cases` |
| `Digenic dominant` | 5 | `digenic_dominant` |
| `Digenic recessive` | 4 | `digenic_recessive` |
| `Oligogenic` | 2 | `oligogenic` |
| `Multifactorial`, `Mitochondrial` | 1, 1 | `multifactorial`, `mitochondrial` |
| the other 17 (`autosomal_dominant`, `somatic`, `x_linked`, `digenic_*`, `imprinting`, `y_linked`, ...) | 5,911 | already the key |

About 1,610 values are dirty, so at most that many records change. Forms that would *not* map - `AD`, `AR`,
`XLR`, the HPO wording `Autosomal dominant inheritance` - are not in the issue's list; the dry run below
reports any that exist on the live data, and only then does Part A gain an alias migration for them.

---

## Data

No model changes. What changes is the value inside one evidence blob:

```python
# before (legacy): a label, or a comma separated string of labels
{"mode_of_inheritance": {"value": "Autosomal recessive, X-linked", "note": "per OMIM", "explain": None}}
# after: the sorted list of option keys, note / explain / immutability untouched
{"mode_of_inheritance": {"value": ["autosomal_recessive", "x_linked"], "note": "per OMIM", "explain": None}}
```

Part C adds to chosen `snpdb/models/models.py:Lab` rows' `classification_config` (merged, not overwritten):

```json
{"mode_of_inheritance": {"mandatory": true}}
```

---

## Part A - `classification_normalise_option_values` (public repo)

A one-off command in `classification/management/commands/` (`category = "one-off"`), general over any
select / multi-select key, since legacy labels are not an SA Path peculiarity:

```bash
python3 manage.py classification_normalise_option_values mode_of_inheritance           # dry run: report only
python3 manage.py classification_normalise_option_values mode_of_inheritance --apply
```

1. **Candidates**: the pks of `Classification.objects.filter(lab__external=False, evidence__<key>__value__isnull=False)`,
   ordered by pk. External labs' records are read-only copies of another instance and are left alone.
2. **Per record**, loaded fresh by pk: take the key through `classification.evidence_keys` (so lab / org option
   overrides apply), split the stored value exactly as `process_entry` does and match each part with
   `Classification.match_option`, dedupe and sort (as `process_option_values` does). Move the split out of
   `process_entry` into a small static helper both call, so the command and every normal write can never
   disagree on what the parts are.
   - Normalised value equals the stored one: **clean**, nothing written.
   - Any part matches no option: **left** untouched as a whole, reported with pk and value.
   - Otherwise, with `--apply`, patch it:
     ```python
     classification.patch_value({key: {"value": normalised, "note": blob.get("note"), "explain": blob.get("explain")}},
                                user=admin_bot(), source=SubmissionSource.VARIANT_GRID, save=True)
     ```
     `VARIANT_GRID` can change an API-immutable value (a `FORM` patch would be refused with an `immutable`
     warning) and `merge_from` keeps the existing immutability. `note` and `explain` go in explicitly
     because a non-form patch with no `explain` sets it to None (the same reason
     `Classification.revalidate` passes both). An `immutable` warning in the response is reported as a failure.
   - Then publish with `classification.publish_latest(user)` (existing share level) **only** if, before the
     patch, `Classification.has_outstanding_changes` was False and the record is not withdrawn. A record with a
     curator's unsubmitted edits gets the normalised value in its draft, published when they submit; a withdrawn
     record is normalised but its withdrawn published version is not re-sent anywhere.
3. Each record is its own `transaction.atomic()`; an exception is logged with the pk and the run carries on.
   Re-running is a no-op for everything already clean.
4. **Report** (both modes): `stored value -> normalised` with counts, then totals per lab for clean /
   normalised and published / normalised, not published (unsubmitted edits, withdrawn) / left (unmatched, with
   pks) / failed. It also prints how many MOI-level `ConditionTextMatch` rows (`classification` null,
   `mode_of_inheritance` set) with non-empty `condition_xrefs` sit on a dirty MOI, since publishing moves each
   classification under the MOI node for its normalised value (see Side effects) - expected to be zero at SA
   Path, where the matching UI is off.

One line in `classification/AGENTS.md` Gotchas: legacy records can hold option labels rather than keys, and
this command is the fix.

### Scale

A plain loop in one process. SA Path, the largest internal classification set, has ~13k internal
classifications, ~6.7k with MOI, and ~1.6k to write. The batching rule in `claude/guides/operations.md` is about
millions of variants / alleles held in memory and hours-long cursors; here the pk list is 6.7k integers and each
record is fetched on its own. A patch plus publish (with its post-publish signals) is well under a second, so
the whole run is minutes, and per-record transactions already give "one failure only loses its own record".
The value tally is a Python counter over those rows because the mapping is `match_option`, which SQL cannot do;
its size is bounded by the number of distinct stored values, not by rows.

### Side effects of republishing

- `classification_post_publish_signal` recalculates the record's `ClinicalContext` - clinical significance
  is unchanged, so no discordance moves.
- `ConditionTextMatch.sync_condition_text_classification` keys the MOI level of the matching hierarchy on the
  stored MOI list (`classification/models/condition_text_matching.py`), so each record moves from the node for
  `Autosomal recessive` to the one for `["autosomal_recessive"]`. A resolution made on the old MOI node would
  not follow; the report counts those. This is why Part A runs before the condition matching backfill.
- Records at a discordant share level get a new last published version, so SA Path's Shariant upload
  (`sync/shariant/variant_grid_upload.py`) re-sends them on its next run - up to ~1.6k, carrying the same
  meaning in the form Shariant already stores.

## Part B - run it on SA Path deployments (sapath repo)

A `ManualOperation` in a new sapath migration (*sapath/migrations/0014_normalise_mode_of_inheritance.py*),
after the pattern of `snpdb/migrations/0188_one_off_migrate_common_filter_gnomad_versions.py`:

```python
ManualOperation(task_id=ManualOperation.task_id_manage(["classification_normalise_option_values", "mode_of_inheritance", "--apply"]),
                note="Convert legacy mode_of_inheritance labels to option keys (SACGF/variantgrid_sapath#363)",
                test=_has_legacy_mode_of_inheritance)
```

`_has_legacy_mode_of_inheritance(apps)` reads the historical `EvidenceKey` option keys and returns True when an
internal-lab classification holds an MOI value (or list part) that is not one of them. The sapath app is only
installed on SA Path deployments, so Shariant and variantgrid.com never register it,
and a fresh or already-clean database does not either.

Rollout on SA Path prod: deploy, and before the migrator's auto-manage pass run the dry run by hand
(read-only). If the report shows no unmatched values and no MOI-level resolutions, let the migrator run the
step; otherwise add the aliases (Decision 1) first.

## Part C - mandatory, after Karin agrees

A stakeholder decision, so nothing here lands until Karin (SA Path) signs off, and it follows Part B on prod.
Then a sapath data migration merges `{"mode_of_inheritance": {"mandatory": true}}` into the
`classification_config` of the germline labs she names, keeping whatever config they already have.

How that behaves (traced, nothing to build):

- `classification/models/classification.py:Classification.evidence_key_overrides` merges org then lab
  `classification_config`, and `classification/models/evidence_key.py:EvidenceKeyMap` sets the override
  attributes on its copy of the key - so only those labs see it, and Shariant (which has its own config for
  the SA Path labs) is unaffected.
- A new record gets a `mandatory` **error** on the empty field: `Classification.create` patches with
  `clear_all_fields=True`, which wipes every mandatory key so `Classification.process_entry` flags it. On the
  form an error replaces Submit with "Resolve Messages to Submit"
  (`variantgrid/static_files/default_static/js/vc_form.js`); through the API
  `BulkClassificationInserter` refuses the share with `share_failure`. So an unfilled record is never shared,
  to the lab or to Shariant.
- An existing record is not touched: validation is stored in each cell when that cell is patched and
  `Classification.has_errors` reads what is stored, so a record already missing MOI only gets the error if
  someone clears MOI or it is revalidated.
- A lab-level override applies to every record in the lab regardless of allele origin (MOI has no namespace),
  which is why only germline labs get it.
- Case report finalising (`classification/report/case_report_builder.py`) publishes without checking errors,
  but it skips records with unsubmitted edits, and a new record without MOI could never have been submitted.

---

## Tests

One new module, *classification/tests/test_classification_normalise_option_values.py*, on
`classification/tests/models/test_utils.py:ClassificationTestUtils` (the MOI key comes from the seed
migration). Legacy records are made by creating and publishing normally, then writing the label straight into
`evidence` and the modification's `published_evidence` with `update()` - the way VG3 left them.

- A comma separated legacy string and a mixed list (`["Isolated cases", "digenic_dominant"]`) become sorted
  option keys, get published, and keep their note and explain.
- The default dry run writes nothing.
- A record with unsubmitted edits is normalised but its last published version does not move.
- A record with one unmatchable part is left entirely as it was.
- An API-immutable value is normalised (the reason for the `VARIANT_GRID` source).
- A clean record gets no new modification.

Nothing tests Part C: it is configuration of existing, tested `EvidenceKeyOverrides` behaviour. Nothing tests
the migration's `test=` function.

## Decisions (defaults - overrule any)

1. **No evidence key change and no aliases**: the list is already defined and every value in the issue maps.
   Aliases only if the SA Path dry run shows forms that don't.
2. **Current evidence only, history untouched.** Old `ClassificationModification`s keep what was submitted then;
   the record gets one new modification by `admin_bot`, which is the audit trail. No note is added to the
   evidence: the change is lossless and a note would travel with the record forever.
3. **`SubmissionSource.VARIANT_GRID`**, not `FORM` as in
   `classification/management/commands/classification_set_legacy_allele_origin.py`, so API-immutable values
   are converted too.
4. **Unsubmitted edits and withdrawn records are normalised but not published.**
5. **A record with any unmatched part is left whole**, not partially converted.
6. **Internal labs only.**
7. **A general command over any select key, dry run by default**, living in the public classification app.
8. **A plain loop**, per-record transaction, carry on after an error.
9. **A `ManualOperation` in the sapath repo**, not in public classification migrations (which would also rewrite
   other organisations' records on Shariant) and not an unrecorded hand run (SA Path has several deployments).
10. **Mandatory per germline lab via a sapath data migration after Karin signs off**, with no sweep that puts
    existing records without MOI into error.
11. **Part A runs before the condition matching backfill.**

## Deferred

- Suggesting an MOI from the condition: the `HP:0000005` descendants related to the condition's term, as traced
  in the issue, or the GenCC gene/disease sources in `OntologyTermRelation.extra` (what
  `ontology/models/models_ontology.py:OntologyVersion.moi_and_submitters` lists for the analysis filter form -
  the whole version's MOI vocabulary, not a per-term lookup, so the helper would be new code). With the field
  mandatory it may not earn its keep.
- A backlog view of existing SA Path records with no MOI (~6.4k at the issue's count), for curators to work
  through.
- Checking Shariant for the same legacy labels; the Part A command would serve if they are there.
