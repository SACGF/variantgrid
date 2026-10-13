# Analysis Status page (rename of Analysis Issues, with scheduler queue graph)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-10-13
Status: in progress

## Goal

The admin page at `/analysis/analysis_issues` ("Analysis Issues") already carries the scheduler's operational brake
(pause / unpause of `snpdb.models.JobsControl`) and a stacked status bar, but that bar counts every node ever made, so
"Finished" dwarfs the rest and the live queue is invisible. Rename the page to **Analysis Status** and add a
scheduler section that answers, at a glance, the questions an operator asks of the node scheduler: how much is
waiting, how much has been dispatched to workers, how much is running, and whether anything is stuck. The graph
style mirrors the per-version stacked bar on the annotation Pipeline Runs page
(`variantgrid/static_files/default_static/js/annotation/variant_annotation_runs.js`).

No database models change.

## Data

A new module `analysis/scheduler_status.py` owns the queries. It exposes one dataclass and two functions.

```python
@dataclass(frozen=True)
class SchedulerStatus:
    waiting: int          # DIRTY, leasable now (no run_after in the future)
    backing_off: int      # DIRTY with a current-version NodeTask whose run_after > now
    dispatched: int       # QUEUED with a live lease - handed to analysis_workers, not yet started
    running: int          # LOADING / LOADING_CACHE with a live lease
    lapsed: int           # QUEUED / LOADING / LOADING_CACHE whose current-version lease is missing or expired
    analyses_with_work: int
    in_flight_target: int  # settings.ANALYSIS_NODE_DISPATCH_BACKLOG_IN_FLIGHT_TARGET
    queues: dict[str, dict]  # {"analysis_workers": {"messages": n, "consumers": n}, "scheduling_single_worker": ...}


@dataclass(frozen=True)
class AnalysisWork:
    analysis_id: int
    analysis_name: str
    username: str
    waiting: int
    dispatched: int
    running: int
    lapsed: int
    oldest_modified: datetime   # oldest `modified` among the analysis's loading nodes
```

Bucket definitions come from the scheduler itself, so the page and the dispatcher agree:
- "Live lease" is the condition `analysis/tasks/analysis_update_tasks.py:_live_lease_q` encodes: the NodeTask on the
  node's *current* NodeVersion has `lease_expires >= now`. A lease on a superseded version says nothing.
- `lapsed` is what `analysis/tasks/analysis_update_tasks.py:lease_ready_nodes` will reclaim (QUEUED) or
  permanently fail (LOADING / LOADING_CACHE) on its next pass - the number an operator most wants to see non-zero
  before a worker is blamed.
- `backing_off` is the `run_after > now` gate in the same function.

### Queries

Aggregation is in SQL. The base queryset is `AnalysisNode.objects.filter(status__in=NodeStatus.LOADING_STATUSES)`
annotated with the current-version NodeTask's `lease_expires` and `run_after` via `Subquery` on
`NodeTask.objects.filter(node_version__node_id=OuterRef("pk"), node_version__version=OuterRef("version"))`. A
subquery annotation rather than a join keeps one row per node, so `Count("pk", filter=Q(...))` needs no `distinct`.
From that one queryset:

- `get_scheduler_status(now)` is a single `.aggregate()` with five conditional counts plus
  `Count("analysis_id", distinct=True)`.
- `get_analyses_with_work(now, limit=50)` is `.values("analysis_id", "analysis__name", "analysis__user__username")`
  with the same conditional counts and `Min("modified")`, ordered by oldest first, sliced to `limit`.

Queue depths reuse `library/vg/status.py`: rename `_queues` to a public `broker_queue_depths(queue_names)` that takes
the names to declare (the `vg status` caller passes every `CELERY_TASK_QUEUES` name; the page passes
`["analysis_workers", "scheduling_single_worker"]`). The passive declare already swallows a broker error into
`{"error": ...}`, so the page renders when RabbitMQ is down - a connection failure on `connection_for_read()` is
caught at the same level and reported the same way for every requested queue.

## Page

### Rename

- View: `view_analysis_issues` becomes `analysis/views/views.py:view_analysis_status`; the template becomes
  `analysis/templates/analysis/view_analysis_status.html`; URL named `analysis_status` at path analysis/analysis_status;
  menu item in `variantgrid/menus.py` titled "Analysis Status". The old path analysis/analysis_issues stays as a
  permanent `RedirectView` to the new name (precedent: `snpdb/urls.py` docs redirect), so bookmarks and the
  `JobsControl.pause` reason text ("Paused from analysis status page by ...") read correctly.
- The failed-node grid keeps its name (`analysis_node_issues_datatable`, `AnalysisNodeIssuesColumns`,
  `AnalysisNodeIssuesFilterForm`): it lists issues, and that is still what it is.

### Layout, top to bottom

1. The pause / unpause alert and buttons, as now, plus a "Dispatch backlog now" button that queues
   `dispatch_analysis_backlog.si().apply_async()` and adds an INFO message (the analogue of the Pipeline Runs page's
   "annotation-scheduler" button). Disabled with a tooltip while jobs are paused, since the task returns immediately.
2. **Scheduler** card: a stacked horizontal Plotly bar in the Pipeline Runs style with, in order,
   Waiting (grey), Backing off (amber), Dispatched (steel blue), Running (green), Lapsed (red) - only live work,
   so the bar is empty when the scheduler is idle and says so in words ("Nothing queued"). Under the bar,
   `labelled` rows for each count, "In flight / target" as `dispatched + running` over `in_flight_target`, and one
   row per broker queue showing messages and consumers (or the error name). The JS moves into a new
   `variantgrid/static_files/default_static/js/analysis/view_analysis_status.js` (named for its template, #2021; the same `showSummaryGraph` shape as
   the annotation one, parameterised by the series list), and the existing all-nodes summary bar is drawn by the
   same function with its current series.
3. **Analyses with outstanding work** table (plain HTML, it is at most `limit` rows): analysis (linked to
   `analysis` view), user, waiting / dispatched / running / lapsed, oldest modified as a timestamp with
   `TableFormat.timestamp`-equivalent `timesince`. Hidden when empty.
4. The existing all-nodes status summary bar and per-status counts, retitled "All nodes".
5. The existing failed-node filter form and grid, unchanged.

### Query budget

`vg page /analysis/analysis_status --queries` before and after. The new section adds three database queries
(scheduler aggregate, per-analysis table, JobsControl) and one broker round trip. Everything else on the page is as
it was.

## Tests

`analysis/tests/test_analysis_status.py`, on `analysis/tests/utils.py:AnalysisSetupMixin`:

- `get_scheduler_status` buckets: one analysis with nodes set directly to each status, NodeTask rows created with
  a live lease, an expired lease, and a future `run_after`, plus a lease on a superseded NodeVersion (which must
  count as lapsed, not running). Assert each bucket and `analyses_with_work`.
- `get_analyses_with_work` returns one row per analysis with the right counts and the oldest `modified`.
- The page renders for a superuser and the old URL redirects to the new one.

Broker depths are patched out in the page test (`broker_queue_depths` returns a fixed dict).

## Docs

`analysis/AGENTS.md` gains one line under "Patterns here" pointing at `analysis/scheduler_status.py` as the
page-facing view of the lease state, with the rule that its buckets are defined in terms of `_live_lease_q` so the
page and the dispatcher cannot disagree. `claude/maps` regenerate via `scripts/vg map` after the URL change.
