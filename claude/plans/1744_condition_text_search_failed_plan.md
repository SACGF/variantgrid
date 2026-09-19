# #1744 — ConditionText: record when the Monarch search failed so automatch can be retried

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

[#1744](https://github.com/SACGF/variantgrid/issues/1744). Sibling of #1742 (make the propagated error
name Monarch), which is partly landed (`98136e191`) and stays separate.

## The problem

When the Monarch search is down during automatch, the condition text is left unmatched - which looks
exactly like Monarch legitimately finding nothing. The curator sees the local OMIM fallback (or nothing)
as the suggestion with no hint that the MONDO search never ran, and nothing ever re-attempts it:
`ConditionTextMatch.sync_all` is only reachable from `classification/management/commands/sync_condition_text_matches.py`
(not in `variantgrid/celery.py` beat), and `classification/models/condition_text_matching.py:ConditionTextMatch.attempt_automatch`
otherwise fires only from `sync_condition_text_classification` when a *new* root or gene level is created
(`published` / `check_for_withdrawn` receivers) or from the admin "Automatch" action
(`classification/admin/condition_text_admin.py`). An outage during an import leaves a batch of texts
permanently unmatched with no way to identify them short of a full `sync_all` sweep, which itself calls
Monarch for every text without a local MONDO match.

### Drift from the issue

- The issue says `attempt_automatch` catches the error. It doesn't see it: the `try/except` around
  `condition_text_search` is in `classification/models/condition_text_matching.py:search_suggestion`,
  which reports it and falls through to `find_local_term(OMIM)`, returning an ordinary
  `ConditionMatchingSuggestion`. `attempt_automatch`'s own `except` only catches other failures. So the
  failure has to be carried out of `search_suggestion` on the suggestion object, not caught in
  `attempt_automatch`.
- The page always re-searches live: `condition_matching_suggestions` calls `top_level_suggestion` whenever
  the root has no terms, so what a curator sees is a fresh search, not the one that failed. The banner
  wording below says that.
- Everything else checked out: `sync_all` has no beat entry; `report_exc_info` sends to Rollbar when
  configured; `ontology/ontology_matching.py` has the same swallow for the interactive picker (with a
  TODO) and no `ConditionText` to attach anything to.

## Data

```python
class ConditionText(TimeStampedModel, GuardianPermissionsMixin):
    normalized_text = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=CASCADE, null=True, blank=True)
    classifications_count = models.IntegerField(default=0)
    classifications_count_outstanding = models.IntegerField(default=0)
    search_failed = models.DateTimeField(null=True, blank=True)
    # when the last automatch's Monarch search raised; None once an automatch completes without one

    class Meta:
        unique_together = ("normalized_text", "lab")
```

New migration classification/migrations/0182_conditiontext_search_failed.py: one `AddField`, no data
migration, no `ManualOperation`. Existing rows start as `None` - a text that was left unmatched by a
past outage is not distinguishable now and stays a job for `sync_condition_text_matches`.

Only the fact and the time are stored, per the issue: the response body decays in minutes and Rollbar
has it; `ConditionText` is a deduped derived row that `sync_all` deletes at `classifications_count == 0`,
not an audit log.

## Code changes

### `classification/models/condition_text_matching.py`

- `ConditionMatchingSuggestion.__init__`: add `self.search_failed = False` ("the Monarch search raised
  while building this; terms, if any, are the local OMIM fallback"). Not included in `as_json` - the page
  reads the stored timestamp, not the live flag.
- `search_suggestion`: in the existing `except Exception` branch keep `report_exc_info()` and remember the
  failure; set `search_failed = True` on whichever suggestion is returned afterwards (the OMIM fallback or
  the empty one). The early `find_local_term(MONDO)` return and `embedded_ids_check` never call Monarch and
  return with the flag `False`.
- `ConditionTextMatch.attempt_automatch`: after `match = top_level_suggestion(...)`, set
  `condition_text.search_failed = now() if match.search_failed else None`. The existing
  `condition_text.save()` at the end of the `try` persists it. Assignment logic is unchanged - a flagged
  match still goes through `is_auto_assignable` (an OMIM fallback can legitimately assign). The outer
  `except Exception` stays as it is and does not touch the field: the field means "the Monarch search
  failed", which is what the retry re-runs.
- `condition_matching_suggestions` is untouched: it runs on GET and must not write. The flag it gets
  back from `top_level_suggestion` is simply ignored there.
- New `ConditionTextMatch.retry_failed_automatch() -> int`: for each `ConditionText` with
  `search_failed__isnull=False`, ordered by `search_failed`, call `attempt_automatch`; return how many
  were cleared. Stop at the first row whose `search_failed` is still set afterwards - Monarch is still
  down, and each attempt costs up to the 60s timeout plus one retry
  (`classification/models/condition_text_search.py`), so one probe per run is the right amount of
  hammering. The remaining rows keep their older timestamp and are first in line next run.

### classification/tasks/condition_text_tasks.py (new file)

```python
@celery.shared_task(queue='db_workers')
def retry_failed_condition_text_automatch():
    cleared = ConditionTextMatch.retry_failed_automatch()
    logging.info("Condition text automatch retry: cleared %d", cleared)
```

Module docstring: owns the beat entry point for #1744; the logic lives on `ConditionTextMatch`.

### `variantgrid/celery.py`

`app.conf.beat_schedule['retry-failed-condition-text-automatch']`, `HOUR_SECS` (raw seconds, per the
timezone note in that file), next to `reconcile-pending-extractions` which is the same shape of
"re-attempt what an outage parked".

### `classification/management/commands/sync_condition_text_matches.py`

Add `--failed`: prints and calls `ConditionTextMatch.retry_failed_automatch()` then returns, before the
full sync. An operator who knows Monarch is back gets the retry now rather than within the hour, without
the full sweep.

### `classification/templates/classification/condition_matching.html`

Below the `Outstanding Classification Records` row, inside the `mt-3` block:

```django
{% if condition_text.search_failed %}
<div class="alert alert-warning mt-2">
    Automatic matching did not complete: the MONDO text search failed
    {{ condition_text.search_failed|timesince }} ago, so nothing was auto-assigned. It is retried hourly;
    the quick suggestions below come from a fresh search and can be applied now.
</div>
{% endif %}
```

`condition_matching_view` already passes `condition_text`; no view change.

### `classification/admin/condition_text_admin.py`

`search_failed` on `list_display`, and a second `ConditionTextStatusFilter` lookup
`("search_failed", "Monarch search failed")` filtering `search_failed__isnull=False`. This is the "which
ones" list from the issue for an operator, without a datatable column.

### `classification/AGENTS.md`

One gotcha line next to the condition resolution note: `search_suggestion` swallows Monarch failures and
falls back to local OMIM; the only trace is `ConditionMatchingSuggestion.search_failed`, which
`attempt_automatch` writes to `ConditionText.search_failed` and the hourly retry task clears.

Then `scripts/vg map` (new field, task, beat entry).

## Tests

New file classification/tests/models/test_condition_text_matching.py, a `TestCase` using
`classification/tests/models/test_utils.py:ClassificationTestUtils.setUp` for the lab, a
`ConditionText` with a root `ConditionTextMatch` (`gene_symbol=None`, `classification=None`) and a text
with no local match (e.g. `"zzz no such syndrome"`; `find_local_term` runs against the seeded test
ontology). Patch `classification.models.condition_text_matching.condition_text_search` - that name, not
the one in `condition_text_search.py`, since the module imports the function.

1. `attempt_automatch` with the search raising `HTTPError`: `search_failed` is set, no terms assigned, no
   exception escapes.
2. Same row, search then returns `[]`: `search_failed` is cleared.
3. `condition_matching_suggestions` with the search raising: `search_failed` stays `None` (the GET rule).
4. `retry_failed_automatch` with two failed rows and the search still raising: the mock is called once and
   both rows keep a timestamp; with the search returning `[]`: both cleared and the return value is 2.

Test 1 and 3 cover the two rules that are easy to break later (surfacing the failure past the swallow;
never writing from the page). Test 4 covers the stop-on-first-failure branch. The migration and the
banner are checked by hand.

## Manual verification

On this box, Monarch reachable:

1. `python3 manage.py migrate classification` (ask first, per AGENTS.md), then
   `python3 manage.py vg page /classification/condition_matching/<pk> --queries` before and after the
   template change - same query count.
2. In `manage.py shell`, pick a `ConditionText`, set `search_failed=now()` and save; load the page: banner
   shows, suggestions still render. Run `manage.py sync_condition_text_matches --failed`; reload: banner
   gone.
3. Simulate the outage: in the shell, patch `condition_text_search` to raise, call
   `ConditionTextMatch.attempt_automatch(ct)`; `ct.search_failed` is set and Rollbar/log has the report.
   Unpatch, run `--failed`; cleared.
4. After a `celeryd_beat` restart (ask first) `celeryd_beat.log` shows the hourly entry scheduled;
   `vg status` shows `db_workers` consuming it.

## Decisions made

- **Surface the failure on the suggestion, not by re-raising.** Re-raising from `search_suggestion` would
  break the page GET and `condition_match_test_view.py`, both of which want the OMIM fallback. A boolean on
  `ConditionMatchingSuggestion` costs one attribute and leaves every caller's behaviour as it is.
- **Cleared whenever an automatch completes without a Monarch failure**, including when the search was
  not needed (local MONDO hit or embedded id). "Cleared on the next successful search" would leave the flag
  set on a text that has since been fully matched by other means.
- **Only the Monarch failure sets it.** `attempt_automatch`'s outer `except` (anything else going wrong)
  keeps its current behaviour; widening the field to "automatch failed for any reason" makes the hourly
  retry re-run things a retry can't fix.
- **One probe per hourly run while Monarch is down** rather than walking every failed row: a sweep during
  an outage is N × up-to-2-minutes for nothing.
- **Both a beat task and a `--failed` command flag**, sharing one method: the task closes "nothing
  re-attempts it"; the flag is four lines and gives an operator the immediate retry.
- **No datatable column or dashboard count.** The admin filter lists the rows; a column on
  `condition_matchings.html` can follow if curators ask for it.

## Definition of done

- `scripts/vg tests --explain` names `classification.tests` (migration) and it passes with `--keepdb`.
- New tasks module has a docstring; `classification/AGENTS.md` gotcha line added.
- `scripts/vg map` refreshed; `scripts/vg docs check` passes.
- This plan's `Status:` updated.
