# Analysis Benchmarking Plan

Admin-only feature for capturing reload-time snapshots of an analysis, storing them with code/environment provenance, and auto-running them on a Celery Beat schedule so we can spot regressions and confirm tuning wins (e.g. `random_page_cost`, #1547 max_af gate, future #1551 sub-cohort caches) without ad-hoc `profile_analysis_nodes` runs.

Motivating context: the `random_page_cost = 1.1` win that closed out #1546 was discovered manually with `profile_analysis_nodes --planner-diagnostic`. The same kind of regression in the other direction — say, a query plan that accidentally falls back to a sequential scan — would today only surface as a user complaint. This feature catches those before deploy by snapshotting periodically.

## Design summary (locked decisions)

- **Two timing values per snapshot, never one.** `wall_seconds` (start of reload to last node done) and `total_node_seconds` (sum of per-node `load_seconds`). The ratio = effective parallelism; deltas where one moves and the other doesn't are diagnostic.
- **Discard-first warmup.** Each benchmark run does a configurable `prewarm_passes` number of throwaway reloads, then the measured one. Warm cache + warm Q-cache state is the user-perceived workload; cold-start variance dwarfs real signal. Default `prewarm_passes = 1`.
- **Snapshot history is per-analysis.** No cross-analysis aggregation in the model layer; the comparison view computes deltas at read time.
- **Per-node breakdown stored alongside.** Without it, the only triage tool for a regression would be re-running. With it, the snapshot itself shows which node moved.
- **Celery Beat drives the auto-bucket.** No cron involvement. A boolean field on `Analysis` opts an analysis into a single weekly bucket; the bucket task fires Saturday 03:00 local and queues per-analysis benchmarks on a dedicated low-concurrency queue. Multi-bucket (different schedules per group) is out of scope for v1; revisit if/when needed.
- **Snapshot creation is non-destructive.** Reload happens against the live analysis but uses the same code path users hit on "reload nodes" — no separate fork.

## Schema preconditions (already landed)

`NodeVersion` and `NodeCount` are `TimeStampedModel` (analysis migration `0101_nodecount_created_nodecount_modified_and_more.py`). Two consequences:

1. **The retrospective inspector** (`get_cpu_and_walltime(analysis_id)` — see comment thread on the issue this plan was filed under) can compute wall time from `max(NodeCount.created) - min(NodeCount.created)` for an analysis. NodeCount cascades from NodeVersion which is regenerated per node-version bump, so the timestamps automatically reflect the latest load batch — no manual invalidation, no cross-load contamination. This works on existing data the moment the next reload happens, no benchmarking infrastructure required.
2. **The snapshot model can sanity-check its explicit timer** against `NodeCount.created` deltas at write time. If the explicit timer says wall=200s but NodeCount.created spread says 50s, something has gone wrong with the wrapper — emit a warning rather than persisting a wrong number.

The snapshot's primary `wall_seconds` field still comes from an explicit timer (most accurate, captures scheduling overhead before the first node starts), not from NodeCount timestamps — but the timestamps are an essential second source of truth.

## Schema

`analysis/models/analysis_benchmark.py` (new file):

```python
class AnalysisBenchmarkSnapshot(TimeStampedModel):
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE, related_name="benchmark_snapshots")
    run_by = models.ForeignKey(User, on_delete=SET_NULL, null=True)
    run_source = models.CharField(max_length=1, choices=[
        ("M", "Manual (admin button)"),
        ("S", "Scheduled (Celery Beat auto bucket)"),
        ("C", "Command line"),
    ], default="M")

    wall_seconds = models.FloatField()
    total_node_seconds = models.FloatField()
    node_count = models.IntegerField()
    prewarm_passes = models.IntegerField(default=1)

    git_hash = models.CharField(max_length=40)
    git_branch = models.TextField()
    code_dirty = models.BooleanField(default=False)

    hostname = models.TextField()
    pg_version = models.TextField()
    pg_settings = models.JSONField(default=dict)  # whitelist: random_page_cost, work_mem,
                                                  # shared_buffers, effective_io_concurrency,
                                                  # effective_cache_size, max_parallel_workers_per_gather

    annotation_version_id = models.IntegerField(null=True)
    cgc_version = models.IntegerField(null=True)  # spot CGC rebuilds between runs

    notes = models.TextField(blank=True)
    error = models.TextField(blank=True)  # set if the reload failed; wall/total may be partial

    class Meta:
        indexes = [models.Index(fields=["analysis", "-created"])]


class AnalysisBenchmarkNodeRow(models.Model):
    snapshot = models.ForeignKey(AnalysisBenchmarkSnapshot, on_delete=CASCADE,
                                  related_name="node_rows")
    node = models.ForeignKey(AnalysisNode, on_delete=SET_NULL, null=True)
    node_pk_at_snapshot = models.IntegerField()  # survives node deletion
    node_type = models.CharField(max_length=64)
    load_seconds = models.FloatField()
    row_count = models.IntegerField(null=True)
    error = models.TextField(blank=True)
    explain_execution_ms = models.FloatField(null=True)  # only when deep=True
    sql_hash = models.CharField(max_length=16, null=True)  # md5 of normalised SQL


class AnalysisAutoBenchmarkOptIn(models.Model):
    """ Opts an analysis into the weekly Celery Beat benchmark bucket.
        Separate model rather than a boolean on Analysis so the audit trail
        (when it was added, by whom, last_run_at) lives next to the data
        without bloating Analysis itself. """
    analysis = models.OneToOneField(Analysis, on_delete=CASCADE, related_name="auto_benchmark")
    added_by = models.ForeignKey(User, on_delete=SET_NULL, null=True)
    added_at = models.DateTimeField(auto_now_add=True)
    last_run_snapshot = models.ForeignKey(AnalysisBenchmarkSnapshot, on_delete=SET_NULL,
                                           null=True, related_name="+")
    notes = models.TextField(blank=True)
```

## Reload mechanics

The benchmark task wraps the existing reload pipeline; it does not reimplement it. Sequence:

1. Capture environment provenance (git, PG version, settings whitelist) once at task start.
2. For `i in range(prewarm_passes)`: invalidate cached counts on every node, fire the analysis-update pipeline, block until done.
3. Final measured run: same invalidation, same update, but `time.perf_counter()` brackets it.
4. After the measured run completes, walk `analysis.analysisnode_set.all()` and snapshot each node's `load_seconds`, `count`, type, plus optional EXPLAIN.

`analysis.update_all_nodes(blocking=True)` does not exist today as a single call — the analysis update is Celery-driven with status polling. Implement a synchronous wrapper that submits the reload and polls until every node lands in `NodeStatus.READY` or `NodeStatus.ERROR`. Existing patterns:

- `AnalysisNode.NodeStatus` enum already has the terminal states.
- `analysis.tasks.node_update_tasks` has the per-node update task.
- The synchronous wrapper just blocks on `apply_async().get()` for each node, or polls `NodeStatus` rows until quiescent.

For the wall-time measurement to be meaningful, the synchronous wrapper has to actually wait — partial completions or returning early would silently corrupt the metric.

## Wall vs total — what they distinguish

| pattern | meaning |
|---|---|
| `wall` and `total` both up | real cost increase — chase the per-node delta |
| `total` flat, `wall` up | contention (locks, autovacuum, busy server) |
| `total` flat, `wall` down | better parallelism — same work, smarter plan |
| `total` down, `wall` flat | parallelism regressed but each piece got faster |
| `total` and `wall` both down | unambiguous win |

The `random_page_cost = 1.1` change manifests as the third row: same work, smaller wall via better plans. The benchmark catches this only because both metrics are stored.

## Cache-state discipline

Single configurable: `prewarm_passes` on the snapshot. Default 1. Recorded on the snapshot row so historical snapshots taken under different disciplines stay distinguishable.

- `prewarm_passes = 0` measures cold start — useful for "post-deploy first user experience" but high variance.
- `prewarm_passes = 1` (default) measures warm cache, fresh Q-cache — typical post-deploy steady state.
- `prewarm_passes = 2+` measures fully-warm — typical mid-session user.

Compare snapshots only against same-prewarm peers in the UI.

## Admin UI

- Add "Benchmark" button next to "Reload Nodes" on the analysis page, visible only when `request.user.is_staff`.
- POST submits a small form: `prewarm_passes`, `deep` (capture EXPLAIN), `notes`. Defaults: 1, false, empty.
- Submit fires `benchmark_analysis_task.apply_async(args=[analysis.pk, options], queue="benchmark_workers")`.
- Redirect to the snapshot detail page. Polling indicator until the task finishes.
- **Snapshot history view** at `/analysis/{pk}/benchmarks/`: table of snapshots sorted by `-created`, columns wall / total / git_hash / source / notes. Each row links to detail.
- **Snapshot detail view**: timing summary, env provenance, per-node table sorted by `load_seconds desc`, plus a "compare to" picker that loads any other snapshot on the same analysis and computes per-node deltas.
- **Sparkline** of `wall_seconds` over the last 30 days on the analysis page header — at-a-glance drift indicator.

## Celery Beat auto bucket

Add to `variantgrid/celery.py` in the existing `app.conf.beat_schedule` dict:

```python
app.conf.beat_schedule['analysis-benchmark-weekly'] = {
    'task': 'analysis.tasks.benchmark_tasks.run_weekly_benchmark_bucket',
    'schedule': crontab(hour=3, minute=0, day_of_week='sat'),
}
```

Saturday 03:00 was chosen to land outside the existing scheduled work that's visible in the project today: the daily DB activity in the early evening (`notify-server-status` at 19:00), the weekly classification email job (Mon 10:00), the seqauto scans (06:00 and 19:00), and the routine RDS-snapshot windows AWS schedules in the early hours of weekdays. Saturday 03:00 is empty as far as the codebase knows.

The bucket task in `analysis/tasks/benchmark_tasks.py`:

```python
@celery.shared_task(queue="benchmark_workers")
def run_weekly_benchmark_bucket():
    """ Fired by Celery Beat. Walks every Analysis with an opt-in row and queues
        a benchmark task per analysis. Sequential by virtue of queue concurrency=1
        on benchmark_workers — they don't compete with each other or with normal
        analysis traffic. """
    for opt_in in AnalysisAutoBenchmarkOptIn.objects.select_related("analysis"):
        benchmark_analysis_task.apply_async(
            args=[opt_in.analysis_id, {"prewarm_passes": 1, "run_source": "S"}],
            queue="benchmark_workers",
        )


@celery.shared_task(queue="benchmark_workers", bind=True)
def benchmark_analysis_task(self, analysis_pk, options):
    """ Performs one snapshot. Called manually (admin button), by the bucket task,
        or by the management command. """
    ...
```

`benchmark_workers` is a new queue with `concurrency=1` so per-analysis benchmarks don't fight each other for shared buffers / CPU. Add to `variantgrid/settings/components/celery_settings.py`. Workers picked up the same way other dedicated queues already are.

Opting an analysis in: small "Add to weekly benchmark" admin action on the analysis page that creates the `AnalysisAutoBenchmarkOptIn` row. Removing it deletes the row.

## Management command (for ad-hoc / CI)

`analysis/management/commands/benchmark_analyses.py`:

```bash
python3 manage.py benchmark_analyses --analyses 21 50 51 \
    --prewarm-passes 1 --notes "post-rpc11-tuning"
```

Runs the same `benchmark_analysis_task` body inline (no Celery hop) for each analysis. Useful for CI sweeps and one-off "did that change help?" measurements.

## Validation / acceptance

1. **Smoke test on dev cohort 6 analysis (#21)**: queue manually, verify snapshot row created with sane `wall` / `total`, per-node rows match `analysis.analysisnode_set` counts.
2. **Re-run discipline**: queue twice in quick succession with `prewarm_passes = 1`. Wall variance < 20% on a quiescent dev DB. If variance is higher, increase the warmup count or investigate cache-state effects.
3. **Bucket task**: opt one analysis in, manually fire `run_weekly_benchmark_bucket()` from a Django shell, confirm snapshot lands with `run_source = "S"` and the opt-in row's `last_run_snapshot` updates.
4. **Compare view**: take two snapshots on the same analysis (with a SQL or schema change in between), open the compare view, confirm per-node deltas match what you expect.
5. **Failure mode**: simulate a node error during reload (e.g. break a sample link) and confirm the snapshot is still written with the partial data and `error` populated, not silently lost.

## Stretch features — not v1

- **Alerting**: if `wall_seconds > 1.5 * trailing_30_day_average`, fire a notification via the existing `library.log_utils.NotificationBuilder` pipeline. Triggered at the end of `benchmark_analysis_task`.
- **Per-node EXPLAIN auto-capture on regression**: when a node's `load_seconds` exceeds its previous snapshot by > 50%, automatically re-run with EXPLAIN ANALYZE and stash the plan as a JSON blob. Targeted forensics, only triggers on regressions, low storage.
- **Settings sweep mode**: `benchmark_analyses --analyses 21 --pg-setting random_page_cost=1.1 --pg-setting random_page_cost=4.0` runs each analysis once per setting variant. Mirrors the `--planner-diagnostic` flag we added to `profile_analysis_nodes`.
- **Multiple auto-bucket schedules**: a M2M `AnalysisBenchmarkBucket(name, schedule_crontab)` model so different analyses can run on different cadences (e.g. critical sentinels nightly, others weekly). Skipped in v1 because a single weekly cadence covers the regression-detection use case; add when there's a concrete need.
- **CI integration**: nightly job runs sentinel analyses on a test DB, fails the build if `wall_seconds` regressed by > 25 % vs the prior `master` snapshot. Catches code regressions before deploy. Belongs in CI config, not the application.

## Files touched

- `analysis/models/analysis_benchmark.py` — new models (snapshot, node row, opt-in).
- `analysis/migrations/00XX_analysis_benchmarking.py` — three `CreateModel`s plus indexes.
- `analysis/tasks/benchmark_tasks.py` — `benchmark_analysis_task`, `run_weekly_benchmark_bucket`, the synchronous reload-wait wrapper.
- `analysis/views/benchmark_views.py` — admin button POST handler, snapshot history view, snapshot detail view, compare view.
- `analysis/templates/analysis/benchmark/{history,detail,compare}.html` — new templates.
- `analysis/urls.py` — three new URL patterns under `/analysis/{pk}/benchmarks/`.
- `analysis/management/commands/benchmark_analyses.py` — CLI entry point.
- `variantgrid/celery.py` — add `analysis-benchmark-weekly` to `app.conf.beat_schedule`.
- `variantgrid/settings/components/celery_settings.py` — add `benchmark_workers` queue.
- `analysis/admin.py` — register the three models for inspection.

## Open items / things to confirm before implementation

1. **Synchronous reload wrapper**. Verify `AnalysisNode.NodeStatus` polling is the right shape, or whether there's a higher-level "reload all and wait" already that I missed in the survey above.
2. **PG settings whitelist**. The snapshot stores a JSON blob of the planner-relevant settings. Confirm the list (`random_page_cost`, `work_mem`, `shared_buffers`, `effective_io_concurrency`, `effective_cache_size`, `max_parallel_workers_per_gather`) is sufficient, or add `jit`, `max_parallel_workers`, etc. if relevant.
3. **`benchmark_workers` queue concurrency**. `1` keeps benchmarks from interfering with each other; verify it doesn't starve the bucket task on a dev box where one analysis benchmark might take 5 minutes and the bucket has 20 analyses opted in. Acceptable answer: 20 × 5 min = 100 min running on Saturday morning, fine.
4. **Anonymous auto-bucket vs UI-restricted**. Bucket runs as no user (`run_by = NULL`); confirm the snapshot detail view handles that correctly (no broken `run_by.email` dereferences).
5. **CGC version capture**. The snapshot's `cgc_version` field is intended to flag "the underlying data changed between snapshots so wall-time deltas may be data-driven, not code-driven." Confirm the source of truth — probably `analysis.cohort_genotype_collection.cohort_version` or similar.
