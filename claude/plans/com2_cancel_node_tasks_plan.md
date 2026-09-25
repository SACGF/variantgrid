# Cancel running node loads when an analysis is deleted or a node bumps version

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-24
Status: landed 5cdebc84b

[SACGF/variantgrid_com#2](https://github.com/SACGF/variantgrid_com/issues/2): deleting an analysis from the analyses
list while one of its nodes is still loading hangs the request ("screen greys out"), or deadlocks in Postgres.

## Why it happens

Deleting from the grid is `snpdb/views/views_permissions.py:group_permissions_object_delete` → `Analysis.delete()`, a
synchronous Django cascade in the web request. The cascade reaches `NodeCache` and
`analysis/models/nodes/analysis_node.py:post_delete_node_cache` drops that node's VariantCollection partition table
(AccessExclusiveLock). A node load mid-query on that partition (a child reading its parent's cache) holds an
AccessShareLock, so the delete blocks behind the load, or the two deadlock.

A node version bump is the same collision one step removed: `AnalysisNode._save` bumps the version, then
`analysis/tasks/node_update_tasks.py:delete_analysis_old_node_versions` deletes the old NodeVersions → same partition
drop → same wait on the old load, which keeps a worker and the database busy until it finishes and only then notices
it is obsolete via `NodeOutOfDateException`.

The per-node cancel already exists: `analysis/views/views_node.py:node_cancel_load` revokes the celery task, aborts it
and runs `pg_cancel_backend(db_pid)`. `analysis/models/nodes/analysis_node.py:NodeTask` records `celery_task` and
`db_pid` in `claim_for_load` / `set_node_task_and_status`. This plan reuses that for the two new triggers.

## Data

No model changes. `NodeTask.celery_task` and `NodeTask.db_pid` are the handles; both are set while a load is running.
Today `_clear_lease` nulls `celery_task` at the end of a load and leaves `db_pid` - see §3.

## 1. Shared cancel helper

New function in `analysis/models/nodes/node_utils.py`:

```python
def cancel_node_tasks(node_task_qs: QuerySet[NodeTask]) -> int:
    """ Stops the celery task and the database query behind each running NodeTask; returns how many were cancelled """
```

For each task in the queryset with `celery_task` set: `app.control.revoke(celery_task, terminate=True)`,
`AbortableAsyncResult(celery_task).abort()`, and if `db_pid` is set `run_sql("select pg_cancel_backend(%s)", [db_pid])`.
Then null `celery_task` and `db_pid` on those rows so a second call is a no-op.

`pg_cancel_backend` is asynchronous, so after signalling, poll `pg_stat_activity` for the cancelled pids still running
a query (`state = 'active'`) for up to `ANALYSIS_NODE_CANCEL_WAIT_SECONDS` (new setting, default 5, in
`variantgrid/settings/components/default_settings.py`) with a short sleep, so a caller that is about to drop partitions
finds the locks released.

`node_cancel_load` becomes a call to this helper for the node's current-version NodeTask, keeping its status update.

## 2. Analysis delete

Landed as `analysis/signals/signal_handlers.py:analysis_pre_delete` (registered in `analysis/apps.py`) rather than in
`pre_delete_analysis`: models_analysis is imported by analysis_node, where NodeTask lives, so the model module cannot
import the helper. It calls `cancel_node_tasks` on
`NodeTask.objects.filter(node_version__node__analysis=instance, celery_task__isnull=False)` before the template
handling, so the cascade's partition drops run against idle backends.

## 3. Version bump

`delete_analysis_old_node_versions` cancels first, then deletes: the tasks to cancel are those on a NodeVersion that is
no longer its node's current version, with `celery_task` set:

```python
NodeTask.objects.filter(node_version__node__analysis_id=analysis_id, celery_task__isnull=False)
    .exclude(node_version__version=F("node_version__node__version"))
```

That task already runs from `analysis/models/nodes/node_utils.py:update_analysis` after every edit and from celery,
so the revoke happens after the edit's transaction has committed.

## 4. Two fixes the helper depends on

- `analysis/tasks/node_update_tasks.py:_clear_lease` also sets `db_pid=None`. Worker database connections persist
  across tasks, so a stale pid on a finished NodeTask would cancel whatever that worker is running now.
- A cancelled query reaches `update_node_task` as `django.db.utils.OperationalError` whose `__cause__` is
  `psycopg.errors.QueryCanceled` (Django 6 on psycopg 3). Today that lands in the transient `OperationalError` branch
  and `_backoff_node` re-queues the load that was just cancelled. Add an `except OperationalError` guard that checks
  for `QueryCanceled` first and exits quietly the way the `NodeOutOfDateException` branch does (the lease is cleared
  in `finally`; for a version bump the new version is already dispatched, for a delete there is nothing left to run).
  Same treatment in `node_cache_task`.

## 5. Tests

In `analysis/tests/test_scheduler.py` (celery is eager under `AnalysisSetupMixin`; patch `app.control.revoke`,
`AbortableAsyncResult.abort` and `run_sql`):

- `cancel_node_tasks` revokes and cancels the backend for a task with both handles, skips one with neither, and nulls
  both handles afterwards.
- `delete_analysis_old_node_versions` cancels the old version's task and leaves the current version's untouched.
- `pre_delete_analysis` cancels every running task in the analysis.
- `update_node_task` with a load raising `OperationalError` from `QueryCanceled` leaves the node's status alone and
  creates no backoff (`run_after` stays null, status is not DIRTY).

## Docs

One gotcha line in `analysis/AGENTS.md` next to the "never call node.save() inside a celery task" entry: a cancelled
query is `QueryCanceled` inside `OperationalError` and means exit, never back off; and `cancel_node_tasks` is the one
place a running load is stopped.
