# AnnotSV import lane: over-committed dispatch + count-task flood

## Symptom

Enabling AnnotSV on `vgaws` (#720 pipeline, `ANNOTATION_ANNOTSV_ENABLED = True`) put the annotation
dispatcher into a state where the DB-import lane barely runs. Observed 2026-08-26 04:20 UTC:

| AnnotSV run status | count |
|---|---|
| FINISHED | 3750 (nearly all empty runs, `dump_count=0`) |
| ANNOTATION_COMPLETED — annotated, waiting to import | **415** |
| CREATED | 176 |
| ERROR | 55 |

AnnotSV itself is healthy on every failed run: stdout ends `...AnnotSV is done with the analysis`,
`pipeline_stderr` empty, return code 0, `.annotated.tsv` written. Nothing is wrong with the tool, the
bundle, or `get_annotsv_command`. The failures are entirely in dispatch and import.

The last non-empty AnnotSV import to complete was 00:30:58; the backlog had been accumulating for ~4
hours. ERROR is the only status any pipeline type has — AnnotSV is not special, it is just the only
pipeline currently producing thousands of short runs, which is what exposes both bugs.

## The three failure shapes

**1. `FileNotFoundError` opening the annotated TSV** (`bulk_annotsv_tsv_inserter.py:173`) — 32 runs.

The run was dispatched to the import lane twice. The first import succeeded (`annotated_count ==
dump_count`), which fired `annotation_run_complete_signal` → `annotation/signals/annotation_run_cleanup.py`
removed the TSV. The second, still-queued import task then opened a file that no longer existed.

Run 43214 is the clean example — `celery_task_logs` holds both attempts:

```
task d3b398e5  import_start 00:45:50.015  import_end 00:45:50.251        <- succeeded
task cf92d646  import_start 01:24:54.985  import_end 01:24:55.084  error <- FileNotFoundError
```

`upload_end` (00:45:50) is *earlier* than `upload_start` (01:24:55) on the row, which is the fingerprint
of this shape. **The data for these runs is in the DB** — only the status is wrong.

**2. "lease expired (worker lost); exceeded max attempts"** with `upload_attempts=0` — 23 runs.

These reached ANNOTATION_COMPLETED with a valid TSV still on disk, were leased three times without their
import task ever reaching a worker, and `reclaim_stalled_annotation_runs` failed them out via
`_fail_stalled_run`. These are genuinely **un-imported** — their AnnotSV columns are missing.

**3. Marked ERROR despite a successful import** — 3 runs (43266, 41662, 43268).

`_fail_stalled_run` wrote `error_exception` while the import task was still queued; the task later ran and
imported cleanly (`upload_end` set, TSV reclaimed), but nothing clears `error_exception`, and
`AnnotationRun.get_status()` (`annotation/models/models.py:1258`) tests it first:

```python
if self.error_exception:
    status = AnnotationStatus.ERROR
```

so a completed upload can never outrank a stale error.

## Root cause

Two independent bugs that compound.

### A — `_dispatch_for_vav` budgets capacity per-VAV, not globally

`_dispatch_sweep` (`annotation_scheduler_task.py:315`) has the correct global budget, and its docstring
names this exact hazard:

> ONE global worker-capacity budget shared across versions (the worker pool is shared, so a per-version
> budget over-commits - each version would see the full pool minus only its own in-flight, and
> collectively lease far more than the pool can run).

`_dispatch_for_vav` (`:380`) — the path `_trigger_dispatch` fires after **every** run completion — never
got that fix. It still computes:

```python
upload_capacity = upload_slots - _lane_in_flight_qs(vav, now, _IMPORT_RUNNING_STATUSES).count()
```

With three ACTIVE VAVs (GRCh37, GRCh38, T2T-CHM13v2.0) each leasing up to
`ANNOTATION_UPLOAD_WORKER_SLOTS = 4`, in-flight reaches 12 against a pool of 4. Confirmed live — import
in-flight was exactly 4 per VAV — and visible in `celery-1.log` on every sweep:

```
dispatch sweep: vep slots=4 in_flight=8 cap=-4 | import slots=4 in_flight=12 cap=-8
```

Negative capacity means the global sweep dispatches nothing at all; only the per-completion path
dispatches, and it over-commits 3×. The surplus tasks sit queued past
`ANNOTATION_RUN_LEASE_SECONDS = 900`, get reclaimed, and burn `attempt_count`.

The VEP lane is over-committed the same way (`in_flight=8` against 4 slots) — it just tolerates it better,
because a VEP run holds its slot with a live heartbeat rather than sitting in a queue.

### B — count kicks are re-issued every dispatch cycle with no in-flight tracking

`_trigger_counts_for_uncounted_runs` (`:211`) is called from *both* dispatch paths and kicks up to
`COUNT_KICK_BATCH = 500` runs per VAV **per invocation**, tracking nothing about what is already queued.
`_trigger_dispatch` fires two dispatch messages per run completion, and with thousands of sub-second
AnnotSV runs finishing that reached ~250 dispatch invocations/minute.

Measured result: ~3000 `count_annotation_run` tasks/minute across all four db_workers, each ~30 ms and
overwhelmingly a no-op re-count of a run already handed out. All four db_workers pegged; import messages
rarely reached one.

`COUNT_KICK_BATCH`'s own comment anticipated the scenario —

> #1646: max count_annotation_run tasks the dispatcher kicks per cycle, so creating a new annotation
> version (thousands of runs at once) doesn't burst the count queue

— the intent was right, but the bound is on *per-cycle* volume rather than *in-flight* work, so at 250
cycles/minute it bursts anyway.

### How they compound

**A** puts more import tasks in the queue than the pool can run; **B** means the queue never drains; so
every import task blows its 900 s lease → reclaim → `attempt_count` burn → either a duplicate import
(shape 1) or a fail-out (shape 2).

## Fixes

### Fix 1 — one global capacity budget

Collapse the two dispatch paths rather than duplicating the budget logic. `_dispatch_sweep` gains a
`priority_vav` that only reorders the VAV list; the budget stays global.

```python
def dispatch_annotation_runs(variant_annotation_version_id=None):
    ...
    vav = VariantAnnotationVersion.objects.get(pk=...) if variant_annotation_version_id else None
    _dispatch_sweep(priority_vav=vav)
```

`_dispatch_for_vav` is deleted. One budget calculation, one place, and the per-completion kick can no
longer over-commit. `_lease_across_vavs` already takes a VAV list and threads remaining capacity through,
so it needs no change beyond Fix 2.

### Fix 2 — make the lease atomic

`_lease_and_launch_run` (`:493`) does an unconditional `.update()`, so two dispatchers racing both "win"
and both `apply_async`:

```python
def _lease_and_launch_run(annotation_run, worker_id, now) -> bool:
    leased = AnnotationRun.objects.filter(
        pk=annotation_run.pk, status=annotation_run.status, task_id__isnull=True,
    ).filter(Q(lease_expires__isnull=True) | Q(lease_expires__lt=now)).update(
        leased_by=worker_id,
        lease_expires=now + timedelta(seconds=settings.ANNOTATION_RUN_LEASE_SECONDS),
        attempt_count=F("attempt_count") + 1,   # see Fix 4 - this moves out
    )
    if not leased:
        return False
    ...apply_async...
    return True
```

Callers decrement capacity only on `True`. The `status=` term also stops a run that finished between the
queryset read and the lease from being launched at all.

### Fix 3 — the import task re-checks state at execution time

This is what actually kills shape 1. A queued `import_annotation_run` (`annotate_variants.py:299`)
currently imports whatever it finds, however stale. After it takes the `task_id` lock:

```python
if annotation_run.status != AnnotationStatus.ANNOTATION_COMPLETED \
        or not (annotation_run.vcf_annotated_filename
                and os.path.exists(annotation_run.vcf_annotated_filename)):
    logging.info("AnnotationRun %s no longer needs importing (status=%s) - stale queued task %s",
                 annotation_run.pk, annotation_run.status, my_task_id)
    return   # finally: releases the lock, writes no error_exception
```

Same predicate as `is_upload_resumable` (`models.py:1369`). A stale duplicate becomes a silent no-op
instead of an ERROR. Note `import_annotsv_tsv` is a pure `bulk_update` over rows the SV VEP run owns, so
it is idempotent — nothing is lost by skipping it, and nothing would have been corrupted by running it.
The file simply wasn't there.

### Fix 4 — `attempt_count` counts executions, not dispatches

Bumped today in `_lease_and_launch_run`, so a run whose task merely *waited in a queue* burns its retry
budget three times and is failed out with a perfectly good TSV on disk (shape 2). That is not what the
setting describes:

```python
# Lease reclaims (dead-worker re-dispatch) allowed before a run is failed to ERROR.
ANNOTATION_MAX_RUN_ATTEMPTS = 3
```

Move the bump to where each task grabs its lock, in both `annotate_variants` and `import_annotation_run`:

```python
num_modified = AnnotationRun.objects.filter(pk=..., task_id__isnull=True).update(
    task_id=my_task_id, attempt_count=F("attempt_count") + 1)
```

A run that never got a worker is then re-dispatched without penalty, which is correct — nothing went
wrong with it.

Trade-off: nothing then caps re-dispatch of a run whose queue is permanently broken. Add a separate
`dispatch_count` (unbounded, log-only) so that condition stays visible rather than silent.

### Fix 5 — counting becomes a third lane in the same sweep

Counting is an optimisation, not a precondition: `_dispatchable_runs_qs` (`:435`) does not filter on
`count`, and a dispatched uncounted-but-empty run just dumps zero variants and finishes. Its value is
avoiding that pointless dump — which matters enormously here, since most ranges hold no SVs, so most
AnnotSV runs are empty.

Dispatch counts the same way runs are dispatched, in the same sweep, immediately before the run lanes.
The key change is that the lane **leases the runs it hands out**, using the lease fields already on the
row. A leased run drops out of the next cycle's candidate set, so 250 sweeps/minute hand out the same
work zero extra times — the flood becomes structurally impossible rather than rate-limited.

```python
def _dispatch_counts(vav, now):
    # One batch in flight per VAV: a leased run isn't a candidate, so the sweep can't re-hand it out.
    already_out = AnnotationRun.objects.filter(
        annotation_range_lock__version=vav, external=False, status=AnnotationStatus.CREATED,
        count__isnull=True, lease_expires__gte=now).exists()
    if already_out:
        return
    run_ids = list(AnnotationRun.objects.filter(
        annotation_range_lock__version=vav, external=False, status=AnnotationStatus.CREATED,
        task_id__isnull=True, count__isnull=True,
    ).filter(Q(lease_expires__isnull=True) | Q(lease_expires__lt=now))
     .order_by("annotation_range_lock__min_variant_id")
     .values_list("pk", flat=True)[:COUNT_KICK_BATCH])
    if not run_ids:
        return
    token = f"count:{uuid4()}"          # mirrors the existing "dispatch:" worker_id convention
    AnnotationRun.objects.filter(pk__in=run_ids).update(leased_by=token, lease_expires=...)
    count_annotation_runs.apply_async((run_ids, token))
```

and in `_dispatch_sweep`, alongside the reclaim already there:

```python
for vav in vavs:
    reclaim_stalled_annotation_runs(vav, now)
    _dispatch_counts(vav, now)
```

**Batch the task rather than kicking 500 messages.** `count_annotation_run(id)` becomes
`count_annotation_runs(ids, token)` looping the existing body. At ~30 ms each, 500 runs is ~15 s on *one*
db_worker instead of 500 messages monopolising all four — the monopoly is what starved the import lane.
The per-run guarded UPDATE swaps its `lease_expires__isnull=True` term for `leased_by=token` ("the lease
is still mine", the shape `save_if_owner` uses), and the task clears the lease on the way out.

What this gets for free:

- **Dead count worker**: `reclaim_stalled_annotation_runs` already matches CREATED + `leased_by` +
  expired lease, and `_reset_run_for_redispatch` on a CREATED run just clears the lease. No new reclaim
  path.
- **No beat, no cache, no new setting, no migration.** `COUNT_KICK_BATCH` survives, now meaning batch
  size.
- **The count/VEP race disappears.** Today an uncounted run can be dispatched to VEP while its count task
  is mid-flight — precisely what the `lease_expires__isnull=True` guard defends against. A real lease
  turns that race into an ordering guarantee.

Two things to get right:

1. **The count lease must not burn `attempt_count`.** `reclaim_stalled_annotation_runs` calls
   `_fail_stalled_run` once `attempt_count >= ANNOTATION_MAX_RUN_ATTEMPTS`, so a run previously
   VEP-dispatched three times could be failed by a *count* lease expiring. Add an explicit branch: a
   `leased_by` starting with `count:` is cleared, never failed and never reset. Fix 4 mostly removes the
   overlap, but the branch should be explicit rather than implied.
2. **A count lease briefly holds runs out of the VEP lane**, since `is_dispatchable` (`models.py:1295`)
   excludes a live lease. With a ~15 s batch that is noise, but it argues for a shorter window than the
   900 s VEP lease — either a small `ANNOTATION_COUNT_LEASE_SECONDS`, or accept the reclaim clearing it.

### Fix 6 — a stale `error_exception` must not outrank a completed upload

Clear `error_exception` on a successful import in `import_annotation_run`. Fix 3 makes shape 3
unreachable, but this is the correct invariant regardless, and it is what makes the recovery below a
one-liner.

### Fix 7 — give the import lane its own queue

`ANNOTATION_UPLOAD_WORKER_SLOTS = 4` is documented as sized to match the db_workers pool "so imports can
use the full pool and drain fast". But `db_workers` is `CELERY_TASK_DEFAULT_QUEUE` — every unrouted task
in the app lands there. The dispatcher is budgeting against a pool it does not own, which is why a flood
of `count_annotation_run` (also db_workers) could starve it to zero.

Route `import_annotation_run` to a dedicated `annotation_import_workers` queue with its own concurrency,
and budget `ANNOTATION_UPLOAD_WORKER_SLOTS` against that. This alone would have contained the outage even
with bugs A and B present. Needs a deployment change (new worker in the celery init config), so it is
correctly a follow-up rather than part of the first pass.

### Optional — debounce `_trigger_dispatch`

`_trigger_dispatch` (`annotate_variants.py:58`) fires two messages per completion, giving ~250
invocations/minute against a single `scheduling_single_worker`. At ~0.26 s each that is ~65 s of work per
minute — saturated, which delays the very sweep that would launch imports.

But this is the same one-off: the rate is driven by thousands of tiny AnnotSV runs completing at once, and
it self-corrects as the backlog drains. A `cache.add(f"dispatch_pending_{vav_id}", ...)` guard collapsing
a burst into one dispatch is cheap insurance, keeping the existing 3 s-delayed second kick for the
commit-visibility race it was written for. Hygiene, not a fix — leave it out of the first pass.

## Ordering

- **First pass — unblocks the running system:** Fix 1 + Fix 5. Together these should drain the 415-run
  backlog.
- **Second pass — closes the races that turned a slow queue into ERRORs:** Fixes 2, 3, 4, 6.
- **Follow-up:** Fix 7, then the debounce if it still looks warranted once the backlog is gone.

## Recovering the 58 currently-errored runs

State re-checked 2026-08-26 04:25 UTC (the count is still growing while the bugs are live):

- **32 `FileNotFoundError` runs** — all have `upload_end` set and `annotated_count == dump_count`. The
  data is already in the DB. Clearing `error_exception` and saving makes `get_status()` recompute to
  FINISHED. No re-annotation needed.
- **23 lease-expired runs** — TSV *and* dump still on disk. `annotation_run_retry(upload_only=True)`
  (`annotate_variants.py:379`) is the right entry point: it already guards on the file existing, and it
  resets `attempt_count` to 0, which matters here.
- **3 lease-expired runs with `upload_end` set** — imported successfully, TSV since reclaimed. Same
  treatment as the first group.

One management command covers all three: clear `error_exception` where `upload_end` is set, call
`annotation_run_retry(upload_only=True)` where the TSV is on disk. Worth adding as a `ManualOperation`
migration (per CLAUDE.md) if other deployments have enabled AnnotSV.

Housekeeping: `/data/annotation_scratch` is at 47 GB / 11,542 files / 439 `annotsv_*` dirs. Errored runs
deliberately keep their output (#1670) and the 415 pending imports keep theirs. Disk is not tight (EFS)
and it should self-reclaim once the lane drains, but it is worth re-checking afterwards.

## Tests

`annotation/tests/test_annotation_dispatch.py` is the home for most of this. Worth keeping:

- Capacity is global across VAVs — three VAVs with pending import runs and 4 slots lease 4 total, not 12.
  This is bug A, and it is the one test that would have caught the outage.
- `_lease_and_launch_run` returns False and launches nothing when the run is already leased / already has
  a `task_id` / has changed status.
- A stale `import_annotation_run` on a FINISHED run is a no-op: no `error_exception`, no second
  `bulk_update`, lock released.
- `attempt_count` is unchanged by a dispatch that never executes, and incremented once by a task that
  does.
- The count lane leases its batch, and a second sweep in the same window hands out nothing.
- A count lease expiring is reclaimed without burning `attempt_count` and without resetting the run.

Not worth keeping: anything asserting that `count_annotation_runs` computes the right count for a given
queryset (that is `get_variants_qs`, already covered), or that `get_status()` maps fields to enum values.

**This host runs against the production DB — do not run `manage.py test` here.** Use `py_compile` for
syntax checking locally and run the suite elsewhere.
