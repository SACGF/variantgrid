# Issue #1413 — API request log records what each request created or changed

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-13
Status: in progress

## Goal

A row on the API Log page (`eventlog/models.py:ApiRequestLog`) says what was sent and what status came back, but
not whether anything was written. The sequencing pipeline re-posts freely, and a re-post is meant to be a no-op
unless something differs, so the question an operator asks is "did this call create or change data, and which
fields?". The log answers it per request: per model, how many rows were created, updated, unchanged and deleted,
and the names of the fields that differed.

Today this cannot be read off the records afterwards: the seqauto serializers go through `get_or_create` /
`update_or_create` at twenty sites and discard the created flag, and `update_or_create` writes unconditionally, so
every re-post bumps `modified` on a `seqauto/models/models_seqauto.py:SeqAutoRecord` whether or not a field
differed. Recording the outcome at write time is the only reliable source, and comparing before writing makes
`modified` mean "changed" as a side effect.

## Data

Two fields on `ApiRequestLog`:

```python
class ApiRequestLog(models.Model):
    ...
    # {"app_label.model": {"created": n, "updated": n, "unchanged": n, "deleted": n, "fields": [names]}}
    # - fields are the ones found different on an update, sorted, across the request
    changes = models.JSONField(default=dict)
    # created + updated + deleted across all models - the grid's sortable "Changes" column and
    # its "wrote data" filter
    num_changes = models.IntegerField(default=0)
```

Existing rows get `{}` and `0`, meaning "not recorded" - the detail page shows nothing for them rather than
"no changes".

The in-request collector, in a new `library/django_utils/change_tracking.py`:

```python
@dataclass
class ModelChanges:
    created: int = 0
    updated: int = 0
    unchanged: int = 0
    deleted: int = 0
    fields: set[str] = field(default_factory=set)


@dataclass
class ChangeCollector:
    by_model: dict[str, ModelChanges] = field(default_factory=dict)  # key "app_label.model"
```

## How changes are collected

`library/django_utils/change_tracking.py` owns a `contextvars.ContextVar` holding the active `ChangeCollector`
(None outside a tracked request), a context manager `collect_changes()` that installs a fresh one and yields it,
and the helper below. A context variable rather than an attribute on the request, because the bulk serializers
are built with `data=` only and have no request in their context.

- **Created / updated** come from a `post_save` receiver connected once at import in
  `library/django_utils/change_tracking.py` (the module is imported by `eventlog/middleware.py`, so it is loaded
  whenever the middleware is). **Deleted** comes from `post_delete`, connected per model by `track_deletes(*models)`
  (`seqauto/apps.py` passes seqauto's models): a receiver with no sender would make `has_listeners` true for every
  model and switch off Django's fast delete across the database. Each receiver is a no-op unless a collector is active. `post_save` records
  `created` from the signal's flag, otherwise `updated`, adding `update_fields` to the model's field names when
  the save passed them. Models in the `eventlog` app are ignored, so the log never counts itself.
- **Unchanged** cannot come from a signal, because the point is that nothing was saved. It comes from
  `update_or_create_if_changed(model, defaults, **lookup)` in the same module: `get_or_create` with the
  defaults, and when the row existed, compare each default against the instance (a foreign key compares by
  `field.attname`, so a model instance and its pk compare equal), `save(update_fields=changed)` only when
  something differs, and record `unchanged` on the collector when nothing does. It returns
  `(instance, outcome)` with `outcome` one of `"created"`, `"updated"`, `"unchanged"`, so callers that branch on
  creation (the sequencing run signal, `set_as_current_sample_sheet`) keep doing so.
- Saves that bypass signals (`bulk_create`, `QuerySet.update`, `QuerySet.delete` on models with no cascade
  receivers) are recorded by the caller with `collector.record(model, created=n)` etc; the seqauto serializers
  have none today, so this is a documented rule rather than a code change.

`eventlog/middleware.py:IntegrationApiMiddleware.__call__` wraps `get_response` in `collect_changes()` and hands
the collector to `_record`, which serialises it into `changes` (fields sorted, empty models left out) and sums
`num_changes`. The `IntegrationActivity.record` upsert happens after the collector is closed.

### Serializer sites

Every `update_or_create` in `seqauto/serializers/*.py` moves to `update_or_create_if_changed`, same lookup and
defaults, with `_created` unpacking replaced by the outcome where it is used. The `get_or_create` sites keep their
semantics (defaults only apply to a new row - `SequencingRun` passes the whole payload as defaults, and a re-post
must not start overwriting it) through the sibling `get_or_create_recorded(model, defaults, **lookup)`, which
records a found row as unchanged. In short: all twenty sites call one of the two helpers. Both take the model
positional-only, so `SequencerModel`'s `model=` lookup passes through.

`Flagstats.objects.create` and other plain creates need nothing - `post_save` counts them.

## Page

- Grid (`eventlog/grids.py:ApiRequestLogColumns`): a "Changes" column rendering `num_changes` with a tooltip of the
  per-model summary, sortable; a "Wrote data" checkbox beside "Errors only" that filters `num_changes > 0`.
- Detail (`eventlog/templates/api_request_log_detail.html`): when `changes` is non-empty, a small table with
  one row per model - verbose name, created, updated, unchanged, deleted - and the changed field names beneath
  the models that have any. A row with `num_changes == 0` but non-empty `changes` is captioned "Re-post, nothing
  differed" so the common case reads at a glance.
- The sequencing run page's API Log tab is the same grid, so it gains the column for free.

## Tests

`eventlog/tests/test_change_tracking.py`:
- `update_or_create_if_changed`: create, update with one differing field (outcome and field name recorded, and
  `modified` bumped), identical re-call (unchanged, `modified` untouched), foreign key given as an instance versus
  a pk compares equal.
- Signals: a plain `objects.create` and a `delete` inside `collect_changes()` are counted; the same outside it
  record nothing; `eventlog` models are ignored.
- Middleware: a tracked POST that creates a row stores `changes` and `num_changes` on the log row; a replay stores
  unchanged with `num_changes == 0`.

`seqauto/tests/`: one existing sample sheet API test extended with a replay assertion - the second identical POST
leaves `modified` alone and its log row shows unchanged.

## Docs

`eventlog/AGENTS.md` (create if absent) gets the rule: a tracked endpoint that writes with `bulk_create` or
`QuerySet.update` records those rows on the collector itself, since no signal fires. `seqauto/AGENTS.md` notes
that serializers write through `update_or_create_if_changed` so a replay is a true no-op and `modified` means
changed.
