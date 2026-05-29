# Analysis Node Scheduling Plan — Issue #346

## Problem

When nodes are added to an analysis while other nodes are already running, the scheduler inserts `wait_for_node` tasks into Celery chains to hold children back until parents finish. The original implementation used `sleep()` inside these tasks, **occupying a Celery worker** for the entire wait. If enough accumulated, the parent nodes could never get a worker — a deadlock.

The VG.com crash log showed exactly this:

```
wait_for_node   Waiting on parent node 78301 which is QUEUED - waiting for 5 secs...
wait_for_node   Waiting on parent node 78300 which is QUEUED - waiting for 5 secs...
wait_for_node   Waiting on parent node 78300 which is QUEUED - waiting for 5 secs...
```

**Root cause** (from the issue): Celery primitives (chords/chains) are fine for a *static* DAG, but an analysis is a *dynamic* DAG — users add/edit nodes mid-run. Embedding blocking `wait_for_node` tasks to bridge a dynamic dependency is fundamentally wrong. As @davmlaw put it:

> *"each worker should only execute 1 task, and then come back, update DB state w/locks, then the next celery worker only takes jobs ready to go… I need to move to determining dependencies and allocating workers using DB locks, so I never assign workers to waiting nodes."*

There is also a second failure mode in the issue body that the previous (purely event-driven) draft of this plan did **not** address:

> *"if the node never finishes (eg celery job was killed and not restarted)… then nothing will ever get done."*

A killed worker leaves the node locked forever, the analysis stalls, and (per another comment) the children are never marked ERROR so the UI looks like it is still working.

---

## Stop-gap Already Applied (variantgrid_com#90)

`wait_for_node` and `wait_for_cache_task` were converted from `sleep()`-based polling to `self.retry(countdown=N)`, releasing the worker between checks. This removed the immediate crash vector but kept the architecture (tasks with unmet dependencies still occupy queue slots, and a killed worker still stalls forever). The stop-gap stays in place; this plan removes these tasks from scheduling chains entirely.

---

## Reference implementation: `mocha/actions`

`mocha` already implements a **state-driven leased job scheduler** backed by the Django DB — the exact shape the issue asks for. This plan ports its design. Key pieces of mocha:

- **`lease_action(worker_id, lease_seconds)`** — claims one ready unit of work using `select_for_update(skip_locked=True)`, sets it to `LEASED`, and stamps `lease_expires = now + lease_seconds`. The ready query is:
  ```python
  Q(status=PENDING) | Q(status=LEASED, lease_expires__lt=now)   # reclaim dead workers
  ```
  plus a `run_after` gate for backoff.
- **`dependencies_satisfied(action)`** — `not action.depends_on.exclude(prerequisite__status=SUCCEEDED).exists()`. A unit is only leasable once every prerequisite has succeeded — *never* run a task with unmet dependencies.
- **`dispatch_pending_actions()`** — the "allocate workers" task: leases every ready action centrally in a loop, then `run_single_action.delay(id)` **one task per leased unit**, fanning the work out across workers in parallel. Doing the leasing centrally avoids every worker hammering the same DB lock.
- **`run_single_action(id)`** — a worker runs exactly **one** pre-leased unit, then returns. On success/permanent-fail it sets terminal state; on transient error it bumps `attempt_count`, applies exponential `next_backoff`, and re-queues via `run_after`.
- **Two trigger sources**: an event path (poll task calls `dispatch_pending_actions.delay()` right after enqueuing work) **and** a periodic loop. The periodic loop is what makes `lease_expires` actually fire — it re-leases units abandoned by dead workers.

### mocha → VariantGrid concept mapping

| mocha | VariantGrid |
|---|---|
| `Action` | `AnalysisNode` (work) + `NodeTask` (lease/lock record) |
| `ActionStatus.PENDING` | `NodeStatus.DIRTY` (waiting to be scheduled) |
| `ActionStatus.LEASED` | `NodeStatus.QUEUED` (claimed by a worker) |
| `ActionStatus.SUCCEEDED` | `NodeStatus.READY` |
| `ActionStatus.FAILED` | `NodeStatus.ERROR` / `ERROR_CONFIGURATION` / `ERROR_WITH_PARENT` / `CANCELLED` |
| `ActionDependency` | `AnalysisEdge` (parent → child) |
| `dependencies_satisfied()` | `_node_ready_to_lease()` — all parents in `READY_STATUSES` **and** required `NodeCache` finished |
| `lease_action()` (`select_for_update(skip_locked=True)`) | `lease_ready_nodes()` |
| `dispatch_pending_actions()` | rewritten `create_and_launch_analysis_tasks()` |
| `run_single_action()` | rewritten `update_node_task()` (one node, then re-trigger) |
| `lease_expires` reclamation | `NodeTask.lease_expires` + periodic sweep |
| periodic `poll` loop | `reschedule_stalled_analyses` Celery beat task |
| `attempt_count` / `next_backoff` / `run_after` | `NodeTask.attempt_count` / `lease_expires` reclaim with backoff |

**The issue's open question — "Make a new state (waiting to be scheduled)? Or is dirty already ok?"** — is answered by mocha: a single waiting state (`DIRTY` ≈ `PENDING`) plus a dependency gate is enough. `QUEUED` becomes the "leased" state. **No new `NodeStatus` value is needed.**

---

## Solution: state-driven leased scheduler

**Core principle (mocha):** a worker only ever runs **one node** for which all dependencies are already satisfied. It never blocks waiting for a parent. After it finishes it updates DB state and re-triggers the dispatcher, which leases whatever just became ready.

```
                 ┌─────────────────────────────────────────────────────┐
                 │  create_and_launch_analysis_tasks(analysis_id)        │
   trigger ─────►│  (scheduling_single_worker)                           │
  (3 sources)    │   lease_ready_nodes(): for each DIRTY/expired node    │
                 │     with all parents READY + caches finished:         │
                 │       claim via select_for_update(skip_locked)        │
                 │       status DIRTY → QUEUED, stamp lease_expires       │
                 │   fan out: update_node_task.delay(node, version) each  │
                 └───────────────┬─────────────────────────────────────┬─┘
                                 │ one task per leased node            │ (nothing ready → fast exit)
                                 ▼
                 ┌─────────────────────────────────────────────────────┐
                 │  update_node_task(node, version)  (analysis_workers)  │
                 │   renew lease; LOADING; node.load()                   │
                 │   → READY / ERROR / CANCELLED; clear lease            │
                 │   finally: trigger create_and_launch_analysis_tasks   │
                 └─────────────────────────────────────────────────────┘
```

**Three trigger sources** (mirroring mocha's event + periodic design, plus VG's existing entry point):

1. **Entry point** — `update_analysis()` after a user edits the analysis (existing).
2. **Event-driven** — every `update_node_task` / `node_cache_task` re-triggers the dispatcher on completion, so a node finishing immediately unblocks its children. *This is the low-latency "kick off allocating workers" behaviour VG wants on top of mocha's loop.*
3. **Periodic safety-net** — a Celery beat task reclaims expired leases (dead workers) and re-dispatches. *This is mocha's periodic loop; it is what makes lease-expiry self-heal.*

### Single-worker invariant (locking/assignment)

**All node leasing/claiming/assignment happens in exactly one place — `lease_ready_nodes`, called only from `create_and_launch_analysis_tasks`, which is routed to `scheduling_single_worker`** (already configured: `CELERY_TASK_ROUTES['analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks'] = SCHEDULING_SINGLE_WORKER`). That queue has a single worker specifically *"to schedule tasks and avoid race conditions"*. Reclaiming abandoned leases and the terminal-fail decision also live in `lease_ready_nodes`, so every state-assignment decision is serialised through that one worker.

Two kinds of DB write are deliberately distinguished:
- **Assignment** (claiming a node for a worker, reclaiming an expired lease, giving up after max attempts) → only ever in `lease_ready_nodes` (single worker).
- **A worker writing the outcome of the node it already owns** (`READY` / `ERROR` / transient-error backoff) → in `update_node_task` / `node_cache_task` on `analysis_workers`, exactly as mocha's `run_action` writes its own action's result. This is not assignment and does not contend on the lease.

Everything else (the entry point, the event re-trigger, the periodic sweep) does no assignment — it only `apply_async`s `create_and_launch_analysis_tasks` into the single worker.

**Why concurrent dispatch is safe** — even though only one worker should normally run the dispatcher, leasing still uses `select_for_update(skip_locked=True)` on a dedicated lock row, exactly like mocha, as belt-and-braces. Two dispatcher runs racing for the same node: one acquires the row lock and transitions `DIRTY → QUEUED`; the other skips the locked row. Combined with `scheduling_single_worker` (single concurrency), triggers for the same analysis serialise with negligible overhead, and the fast-exit path returns when nothing is ready.

---

## Files Changed

| File | Change |
|------|--------|
| `analysis/models/nodes/analysis_node.py` | Add lease/backoff fields to `NodeTask` (`lease_expires`, `leased_by`, `attempt_count`, `last_attempt`, `run_after`) |
| `analysis/migrations/XXXX_nodetask_lease.py` | New migration for the `NodeTask` lease/backoff fields |
| `analysis/tasks/analysis_update_tasks.py` | Rewrite `create_and_launch_analysis_tasks` as lease-and-fan-out; add `lease_ready_nodes` + `_node_ready_to_lease` (reclaim + terminal-fail live here, in the single worker); drop `wait_for_node`, `_get_node_dependencies`, `_get_ancestor_set`, `dependencies` param, group/chord building |
| `analysis/tasks/node_update_tasks.py` | `update_node_task` runs one node, writes its own outcome (incl. transient-error backoff), clears lease, re-triggers; `node_cache_task` re-triggers; add `next_backoff` + `reschedule_stalled_analyses` beat task; eliminate `WAIT_FOR_CACHE_TASK` from scheduling |
| `analysis/models/nodes/analysis_node.py` (`get_cache_task_args_set`) | Treat an in-flight cache as an unmet dependency rather than emitting `WAIT_FOR_CACHE_TASK` |
| `variantgrid/celery.py` | Register `reschedule_stalled_analyses` in `app.conf.beat_schedule` |
| `variantgrid/settings/components/celery_settings.py` | Route `reschedule_stalled_analyses` to `DB_WORKERS` (it only discovers + kicks the single-worker dispatcher) |

---

## Detailed Changes

### 1. `NodeTask` — add lease fields

`NodeTask` already exists as the per-`(node, version)` lock (`unique_together = ("node", "version")`, with `celery_task` and `db_pid`). Extend it into mocha's lease record:

```python
# analysis/models/nodes/analysis_node.py
class NodeTask(TimeStampedModel):
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    version = models.IntegerField(null=False)
    analysis_update_uuid = models.UUIDField()
    celery_task = models.CharField(max_length=36, null=True)
    db_pid = models.IntegerField(null=True)
    # NEW (mocha lease + backoff fields):
    leased_by = models.CharField(max_length=64, null=True)
    lease_expires = models.DateTimeField(null=True)
    attempt_count = models.IntegerField(default=0)
    last_attempt = models.DateTimeField(null=True)
    run_after = models.DateTimeField(null=True)  # backoff gate: not leasable before this

    class Meta:
        unique_together = ("node", "version")
```

Migration `analysis/migrations/XXXX_nodetask_lease.py` adds the five fields (`attempt_count` defaults to 0, the rest nullable — no data backfill needed).

> The `unique_together("node", "version")` row **is** the lease, exactly as mocha's `Action` row is. We keep the existing insert-to-claim semantics for brand-new work and add lease-expiry reclamation for abandoned work.

---

### 2. `_node_ready_to_lease` — the dependency gate (mocha `dependencies_satisfied`)

A node may be leased only when every dependency is satisfied: all parents settled, and any `NodeCache` it needs has finished building.

```python
# analysis/tasks/analysis_update_tasks.py

def _node_ready_to_lease(node, parents_by_child, nodes_by_id) -> bool:
    """ mocha dependencies_satisfied() equivalent.
        All parents must be in READY_STATUSES (settled — READY or any ERROR),
        i.e. none still loading. ERROR parents count as 'ready' so the child
        runs and fails fast with ERROR_WITH_PARENT (matches issue comment). """
    for parent_id in parents_by_child.get(node.id, ()):
        parent = nodes_by_id.get(parent_id)
        if parent is None or NodeStatus.is_loading(parent.status):
            return False
    return True
```

`NodeStatus.LOADING_STATUSES = [DIRTY, QUEUED, LOADING_CACHE, LOADING]`; `READY_STATUSES` includes `READY` and all error states. So a parent is a blocker iff it is still loading.

The **cache** dependency (replacing `WAIT_FOR_CACHE_TASK`) is handled in `get_cache_task_args_set` — see §6.

---

### 3. `lease_ready_nodes` — claim work atomically (mocha `lease_action`)

**Runs only in `scheduling_single_worker`** (it is called only from `create_and_launch_analysis_tasks`). Because the graph is small (tens of nodes) and it runs single-threaded, the SQL filter is kept coarse and the readiness/backoff/reclaim decisions are refined in Python — clearer than expressing parent/cache/run_after gating in one query.

```python
# analysis/tasks/analysis_update_tasks.py

LEASE_SECONDS = MINUTE_SECS * 10   # a node load that exceeds this is already a problem;
                                   # reclaiming it (risking duplicate work) is acceptable. No heartbeat.

def lease_ready_nodes(analysis_id, worker_id, lease_seconds=LEASE_SECONDS):
    """ mocha lease_action() + dispatch, batched per analysis. SINGLE-WORKER ONLY.
        Claims fresh work, reclaims expired leases, and gives up (terminal fail) past
        MAX_NODE_ATTEMPTS. Returns the (node_id, version) tuples now owned by this run. """
    now = timezone.now()
    leased = []
    with transaction.atomic():
        # Coarse base-table filter (no select_subclasses — keeps select_for_update simple):
        # candidates are nodes that are either waiting (DIRTY) or look stalled (loading).
        # skip_locked => a racing dispatcher simply skips rows we hold.
        candidates = list(
            AnalysisNode.objects
            .select_for_update(skip_locked=True)
            .filter(analysis_id=analysis_id, status__in=NodeStatus.LOADING_STATUSES)
            .order_by("pk")
        )
        nodes_by_id, parents_by_child = _load_graph(analysis_id)   # status snapshot + edges
        tasks_by_node = {}                                         # current-version lease rows only
        for nt in NodeTask.objects.filter(node__analysis_id=analysis_id):
            node = nodes_by_id.get(nt.node_id)
            if node and node.version == nt.version:
                tasks_by_node[nt.node_id] = nt

        for node in candidates:
            node_task = tasks_by_node.get(node.pk)
            lease_live = node_task and node_task.lease_expires and node_task.lease_expires >= now
            if node.status != NodeStatus.DIRTY and lease_live:
                continue  # someone is actively working it, lease still valid → leave alone

            # Reclaim path: a non-DIRTY node with an expired (or missing) lease was abandoned.
            attempts = node_task.attempt_count if node_task else 0
            if attempts >= MAX_NODE_ATTEMPTS:
                _fail_node(node, "Node task did not complete (worker lost); marked failed.")
                continue
            if node_task and node_task.run_after and node_task.run_after > now:
                continue  # backing off — not leasable yet
            if not _node_ready_to_lease(node, parents_by_child, nodes_by_id):
                continue  # parents still loading / required cache still building

            NodeTask.objects.update_or_create(
                node_id=node.pk, version=node.version,
                defaults={
                    "analysis_update_uuid": uuid.uuid4(),
                    "leased_by": worker_id,
                    "lease_expires": now + timedelta(seconds=lease_seconds),
                    "last_attempt": now,
                    "attempt_count": (attempts + 1),
                    "run_after": None,
                    "celery_task": None,
                },
            )
            AnalysisNode.objects.filter(pk=node.pk, version=node.version).update(status=NodeStatus.QUEUED)
            leased.append((node.pk, node.version))
    return leased


def _fail_node(node, message):
    with disable_auditlog():
        AnalysisNode.objects.filter(pk=node.pk, version=node.version).update(
            status=NodeStatus.ERROR, shadow_color=NodeColors.ERROR, errors=[{"message": message}])
    NodeTask.objects.filter(node_id=node.pk, version=node.version).update(lease_expires=None, leased_by=None)
```

Notes:
- `status__in=LOADING_STATUSES` (= `DIRTY, QUEUED, LOADING_CACHE, LOADING`) is the candidate set: anything not yet settled. A `QUEUED`/`LOADING` node with a *live* lease is skipped (its worker is busy); one with an expired/missing lease is reclaimed.
- **Reclaim + terminal-fail are here, in the single worker** — not in the beat task — so all assignment is serialised. A node that has burned `MAX_NODE_ATTEMPTS` is failed (`_fail_node`), which makes its children leasable and they settle to `ERROR_WITH_PARENT`.
- `run_after` gates backoff (set by the transient-error path in §5): a DIRTY node whose `run_after` is still in the future is left for a later pass.
- `attempt_count` is incremented on each lease (each "we tried to run it"), bounding both transient-error retries and dead-worker reclaims with one counter.

`_load_graph` reuses `get_nodes_by_id` + the edge query to build `nodes_by_id` (current statuses) and `parents_by_child` (the `parent_value_data` mapping currently built inline in `_get_analysis_update_tasks`).

---

### 4. `create_and_launch_analysis_tasks` — dispatch by fan-out (mocha `dispatch_pending_actions`)

The whole `_get_analysis_update_tasks` / topo-sort / group / chord machinery is replaced by central leasing + one task per node. **No chains, no groups, no `wait_for_node`.**

```python
# analysis/tasks/analysis_update_tasks.py

@celery.shared_task
def create_and_launch_analysis_tasks(analysis_id, run_async=True):
    """ Routed to scheduling_single_worker (the only caller of lease_ready_nodes).
        mocha dispatch_pending_actions(): lease every node ready right now, fan out one task each. """
    worker_id = f"dispatch:{create_and_launch_analysis_tasks.request.id or 'sync'}"
    leased = lease_ready_nodes(analysis_id, worker_id)
    if not leased:
        return  # fast exit — nothing ready (mocha 'nothing to do')

    for node_id, version in leased:
        sig = _node_launch_signature(node_id, version)  # update task, or [cache, update] for cache-building nodes
        if run_async:
            sig.apply_async()
        else:
            result = sig.apply()
            if not result.successful():
                raise Exception(result.result)
```

`_node_launch_signature` is the per-node replacement for `_add_jobs_for_group`: for an ordinary node it returns `node.get_update_task()`; for a node that must build **its own** cache first it returns the two-task chain `chain(node_cache_task.si(...), update_node_task.si(...))` (a single-node chain — no cross-node dependency, so no deadlock risk). Cache builds shared across nodes are deduped by `NodeCache.get_or_create_for_node` + the in-flight-cache dependency gate, not by a wait task.

```python
def _node_launch_signature(node_id, version):
    node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
    update_job = node.get_update_task()
    cache_jobs = [Signature(t, args=a, immutable=True)
                  for (t, a) in node.get_cache_task_args_set() if t]
    if cache_jobs:
        return chain(*cache_jobs, update_job)   # same node: cache then update
    return update_job
```

Each leased node is dispatched independently → multiple `analysis_workers` pick them up in parallel → response times drop. A node finishing re-triggers the dispatcher (§5), which leases the next layer. Topological ordering is no longer pre-computed; it **emerges** from the readiness gate, round by round — which is what lets the DAG change underneath us safely.

---

### 5. `update_node_task` — run one node, then re-trigger (mocha `run_single_action`)

Each worker runs exactly one node, renews then clears its lease, and on every terminal outcome re-triggers the dispatcher so newly-unblocked children get leased.

```python
# analysis/tasks/node_update_tasks.py

@celery.shared_task(base=AbortableTask)
def update_node_task(node_id, version):
    analysis_id = None
    with disable_auditlog():
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
            analysis_id = node.analysis_id            # capture before load()/deletion
        except AnalysisNode.DoesNotExist:
            return                                    # stale / deleted — nothing to re-trigger

        try:
            errors = None
            node_errors = node.get_errors()
            if not node_errors:
                try:
                    node.set_node_task_and_status(update_node_task.request.id, NodeStatus.LOADING)
                    node.load()
                    if node.shadow_color == NodeColors.ERROR and node.is_valid:
                        node.update(shadow_color=None)
                except NodeOutOfDateException:
                    logging.warning("Node %d/%d out of date - exiting", node.pk, node.version)
                    return  # version bumped — reload_analysis_nodes already re-triggered; do NOT re-trigger here
                except OperationalError:
                    # Transient (DB blip / lock timeout). VG node loads are local DB work,
                    # not network calls, so this is rare and usually clears fast — back off
                    # and retry rather than give up, until MAX_NODE_ATTEMPTS.
                    if _backoff_node(node_id, version, analysis_id):
                        return  # node set DIRTY + run_after; delayed re-trigger scheduled
                    status = NodeStatus.CANCELLED  # out of retries
                except NodeConfigurationException:
                    status = NodeStatus.ERROR_CONFIGURATION
                except NodeParentErrorsException:
                    status = NodeStatus.ERROR_WITH_PARENT
                except Exception:
                    errors = get_traceback()
                    status = NodeStatus.ERROR
                else:
                    status = None  # load() set READY itself
            else:
                status = AnalysisNode.get_status_from_errors(node_errors)

            if status is not None:
                try:
                    shadow_color = NodeColors.ERROR if NodeStatus.is_error(status) else None
                    node.update(status=status, errors=errors, shadow_color=shadow_color)
                except (IntegrityError, NodeOutOfDateException):
                    pass  # out of date / deleted
        finally:
            _clear_lease(node_id, version)

    if analysis_id is not None:
        _trigger_rescheduling(analysis_id)


def _clear_lease(node_id, version):
    NodeTask.objects.filter(node_id=node_id, version=version).update(
        lease_expires=None, leased_by=None, celery_task=None)


def _trigger_rescheduling(analysis_id):
    """ Event-driven kick (mocha: poll calls dispatch_pending_actions.delay()).
        Immediate trigger handles the common path; the short-delay trigger covers the
        race where this task commits just after a concurrent dispatcher read statuses.
        Both are no-ops if nothing is ready (lease + fast-exit), and serialise through
        scheduling_single_worker. """
    sig = Signature("analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks",
                    args=(analysis_id,))
    sig.apply_async()
    sig.apply_async(countdown=3)


MAX_NODE_ATTEMPTS = 3

def next_backoff(attempts):
    """ mocha next_backoff(), shorter ceiling. VG node loads are local DB work, not
        network calls, so transient failures are rare and clear fast: 10s, 20s, 40s … cap 5 min. """
    return min(10 * (2 ** attempts), MINUTE_SECS * 5)


def _backoff_node(node_id, version, analysis_id) -> bool:
    """ The worker writes the outcome of its OWN node (not assignment): on a transient
        error, set it back to DIRTY with a future run_after so the single-worker dispatcher
        re-leases it after the delay. Returns False once attempts are exhausted. """
    now = timezone.now()
    node_task = NodeTask.objects.filter(node_id=node_id, version=version).first()
    attempts = node_task.attempt_count if node_task else MAX_NODE_ATTEMPTS
    if attempts >= MAX_NODE_ATTEMPTS:
        return False
    delay = next_backoff(attempts)
    NodeTask.objects.filter(node_id=node_id, version=version).update(
        run_after=now + timedelta(seconds=delay), lease_expires=None, leased_by=None, celery_task=None)
    with disable_auditlog():
        AnalysisNode.objects.filter(pk=node_id, version=version).update(status=NodeStatus.DIRTY)
    # Prompt re-lease after the backoff window (periodic sweep is the catch-all safety net).
    Signature("analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks",
              args=(analysis_id,)).apply_async(countdown=delay)
    return True
```

> `MAX_NODE_ATTEMPTS` / `next_backoff` are shared by both retry reasons: transient-error backoff here, and dead-worker reclaim in `lease_ready_nodes` (§3). One `attempt_count` bounds both. (`MAX_NODE_ATTEMPTS` is referenced from `analysis_update_tasks` too; define it once and import, to avoid duplication.)

The `NodeOutOfDateException` early return is the one path that does **not** re-trigger: the version was bumped by a user edit, which already routed through `reload_analysis_nodes → update_analysis → create_and_launch_analysis_tasks` for the new version.

---

### 6. `node_cache_task` and the end of `WAIT_FOR_CACHE_TASK`

Two changes:

**(a) Treat an in-flight shared cache as an unmet dependency, not a wait task.** Today `get_cache_task_args_set` emits `WAIT_FOR_CACHE_TASK` when another node is already building the cache this node needs. Instead, return *nothing to launch* and let the readiness gate hold the node back; the cache-builder's completion re-triggers the dispatcher.

```python
# analysis/models/nodes/analysis_node.py  get_cache_task_args_set()
    node_cache, created = NodeCache.get_or_create_for_node(self)
    if created:
        task_args_set.add((self.NODE_CACHE_TASK, (self.pk, self.version)))
    # else: cache already being built by another node — do NOT add WAIT_FOR_CACHE_TASK.
    #       _node_ready_to_lease() keeps this node DIRTY until the cache collection
    #       reaches a finished ProcessingStatus; node_cache_task re-triggers on completion.
```

`_node_ready_to_lease` gains a cache check: if `NodeCache.get_or_create_for_node(node)` (without creating) maps to a `variant_collection` whose status is `CREATED`/`PROCESSING`, the node is not yet leasable.

**(b) `node_cache_task` re-triggers on completion** so nodes blocked only on the cache get leased:

```python
# analysis/tasks/node_update_tasks.py
@celery.shared_task
def node_cache_task(node_id, version):
    analysis_id = None
    try:
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
            analysis_id = node.analysis_id
        except AnalysisNode.DoesNotExist:
            return
        # ... existing body unchanged (write_cache, set ProcessingStatus, set node status) ...
    finally:
        if analysis_id is not None:
            _trigger_rescheduling(analysis_id)
```

The single-node `[node_cache_task, update_node_task]` chain (for a node building **its own** cache) is preserved by `_node_launch_signature` — both tasks are the same node, so chaining them is not a cross-node dependency and cannot deadlock.

---

### 7. Lease-expiry reclamation + periodic sweep (mocha's self-healing loop)

This is the part the old plan lacked — it directly fixes *"if the node never finishes (eg celery job was killed)… nothing will ever get done."*

A Celery beat task runs every 60s. It does **no assignment** — it only *discovers* analyses that have schedulable-but-unscheduled work (an expired lease, or a DIRTY node whose backoff `run_after` has passed) and kicks the single-worker dispatcher for each. The actual reclaim / re-lease / terminal-fail all happen inside `lease_ready_nodes` (§3) in the single worker.

```python
# analysis/tasks/node_update_tasks.py

@celery.shared_task
def reschedule_stalled_analyses():
    """ mocha's periodic dispatch loop. Discovery only — finds analyses with work that
        should be runnable but isn't currently leased, and kicks the single-worker
        dispatcher. This is what makes lease-expiry self-heal after a worker is killed. """
    now = timezone.now()
    stalled_lease = Q(nodetask__lease_expires__lt=now)            # abandoned by dead worker
    ready_backoff = Q(status=NodeStatus.DIRTY) & (               # backoff window elapsed (or none)
        Q(nodetask__run_after__isnull=True) | Q(nodetask__run_after__lte=now))
    analysis_ids = (AnalysisNode.objects
                    .filter(status__in=NodeStatus.LOADING_STATUSES)
                    .filter(stalled_lease | ready_backoff)
                    .values_list("analysis_id", flat=True).distinct())

    for analysis_id in analysis_ids:
        Signature("analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks",
                  args=(analysis_id,)).apply_async()  # → scheduling_single_worker (fast-exit if nothing ready)
```

Register the beat task and route the *discovery* task off the single worker (the dispatch it kicks still lands on the single worker):

```python
# variantgrid/celery.py
app.conf.beat_schedule['reschedule_stalled_analyses'] = {
    'task': 'analysis.tasks.node_update_tasks.reschedule_stalled_analyses',
    'schedule': 60,  # seconds; existing comment re: crontab timezone issues — use raw seconds
}
```
```python
# variantgrid/settings/components/celery_settings.py — CELERY_TASK_ROUTES
'analysis.tasks.node_update_tasks.reschedule_stalled_analyses': DB_WORKERS,
```

---

### 8. Removals

```python
# analysis/tasks/analysis_update_tasks.py — DELETE:
from analysis.tasks.node_update_tasks import wait_for_node   # import
def _get_analysis_update_tasks(...): ...                     # replaced by lease_ready_nodes + dispatch
def _add_jobs_for_group(...): ...                            # replaced by _node_launch_signature
def _get_celery_workflow_task(groups): ...                   # no more groups/chains across nodes
def _get_ancestor_set(node_id, parent_value_data): ...
def _get_node_dependencies(nodes_by_id, parent_value_data): ...
```

`wait_for_node` and `wait_for_cache_task` **remain defined** in `node_update_tasks.py` — they are no longer added to scheduling chains, but `wait_for_node` is still called synchronously by `analysis_grid_export_tasks._wait_for_output_node()` as a blocking export gate (a legitimate single-caller use, not a scheduling dependency).

---

## State Transitions (new model)

```
User edits analysis
 └─ reload_analysis_nodes()                 all nodes → DIRTY (version++), NodeVersion rows created
     └─ update_analysis(analysis_id)
         └─ create_and_launch_analysis_tasks (scheduling_single_worker)   ◄── trigger 1: entry
             └─ lease_ready_nodes():  DIRTY nodes with all-parents-READY + caches-finished
                   select_for_update(skip_locked) → status QUEUED, stamp lease_expires
             └─ fan out: update_node_task.delay(node, version)   (one per node, parallel)
                     │
                     ▼
             update_node_task (analysis_workers — writes its OWN node's outcome)
                   renew lease → LOADING → load() → READY / ERROR / CANCELLED → clear lease
                     │                              └─ OperationalError → backoff: DIRTY + run_after,
                     │                                 schedule delayed re-trigger (≤ MAX_NODE_ATTEMPTS)
                     ▼
             _trigger_rescheduling(analysis_id)                  ◄── trigger 2: event-driven
                   create_and_launch_analysis_tasks ×2 (immediate + 3s)  → scheduling_single_worker
                     └─ lease newly-unblocked children → … (repeat until none ready → fast exit)

(independently)  reschedule_stalled_analyses  (Celery beat, 60s, db_workers)  ◄── trigger 3: periodic safety-net
                   discovery only: finds analyses with expired leases / elapsed-backoff DIRTY nodes
                   → kicks create_and_launch_analysis_tasks (single worker), where reclaim / re-lease /
                     terminal-fail (past MAX_NODE_ATTEMPTS → ERROR → children ERROR_WITH_PARENT) happen
```

All claiming/assignment (the `select_for_update` lease, reclaim, terminal-fail) occurs only in `lease_ready_nodes`, reached only via `create_and_launch_analysis_tasks` on `scheduling_single_worker`. Workers on `analysis_workers` only ever write the outcome of the one node they already hold.

---

## What is NOT changed

- **`wait_for_node` definition** — kept for `analysis_grid_export_tasks._wait_for_output_node()` (synchronous export gate). The stop-gap `self.retry()` form is preserved. It is simply never added to scheduling chains.
- **`wait_for_cache_task` definition** — kept (harmless) but never scheduled; in-flight caches are now a dependency gate.
- **`reload_analysis_nodes` / `update_analysis`** — unchanged entry points; they still set nodes DIRTY and call `create_and_launch_analysis_tasks`.
- **`NodeCache` / `node_cache_task` body** — cache *building* logic unchanged; only the trigger-on-completion and the dependency-gate-instead-of-wait change.
- **`scheduling_single_worker`** — still the serialisation point for all dispatch/reclaim.
- **`delete_analysis_old_node_versions`** — unchanged; still cleans stale `NodeVersion` (and cascaded `NodeTask`) rows.

---

## Testing

Acceptance tests already merged to master (per issue): `analysis/tests/test_scheduler.py` has two `@expectedFailure` tests documenting the broken dynamic-dependency behaviour, and `analysis/tests/utils.py` provides `AnalysisSetupMixin`. Removing the `@expectedFailure` markers is the bar for this fix.

1. **`_node_ready_to_lease`** — parent LOADING → False; parent READY → True; parent ERROR → True (child runs, fails `ERROR_WITH_PARENT`); required cache `PROCESSING` → False; cache `SUCCESS`/`SKIPPED` → True.

2. **`lease_ready_nodes` atomicity** — two concurrent dispatcher calls for the same analysis lease each node exactly once (assert no node dispatched twice); the loser leases nothing and the dispatcher fast-exits.

3. **Cascade scheduling** — chain A→B→C, all DIRTY. Dispatch → only A leased/dispatched. Complete A (READY) → re-trigger → B leased, C not. Complete B → C leased. Assert no `wait_for_node` signature is ever produced.

4. **In-flight shared cache** — two MergeNodes sharing a Venn parent's cache: first leases and builds (`node_cache_task`), second stays DIRTY until cache `SUCCESS`, then is leased by the cache task's re-trigger. No `wait_for_cache_task` scheduled.

5. **Dead-worker reclamation** — lease a node, manually set `lease_expires` in the past without completing the task, leave status QUEUED. Run `reschedule_stalled_analyses` (→ dispatcher): node is re-leased (attempt_count increments). After `MAX_NODE_ATTEMPTS`, node → ERROR and a child → ERROR_WITH_PARENT (no infinite hang).

5b. **Transient-error backoff** — make `node.load()` raise `OperationalError` on attempt 1–2, succeed on 3. Assert: after each failure the node is DIRTY with a future `run_after` and `attempt_count` increments; a dispatcher run before `run_after` does **not** re-lease it; after `run_after` it is re-leased and finally reaches READY. With failure on all attempts, it ends CANCELLED after `MAX_NODE_ATTEMPTS`. Assert `next_backoff` produces 10/20/40s.

5c. **Single-worker invariant** — assert `lease_ready_nodes` is invoked only from `create_and_launch_analysis_tasks`, and that `update_node_task` / `node_cache_task` / `reschedule_stalled_analyses` never call it (grep-level or via a guard/mock in tests).

6. **Race — parent finishes just before a dispatch read** — parent commits READY after the dispatcher's status snapshot; assert the 3s-delayed `_trigger_rescheduling` (or the next beat sweep) leases the child.

7. **Stress (manual)** — the auto-run analysis scenario from the crash report and `variantgrid_sapath#175`: all nodes reach READY/ERROR with no worker starvation and no `wait_for_node` accumulation.

---

## Migration / Deployment Notes

- **One DB migration**: `NodeTask` lease/backoff fields (`lease_expires`, `leased_by`, `attempt_count`, `last_attempt`, `run_after`), all nullable / defaulted — no backfill.
- The variantgrid_com#90 stop-gap is already deployed and is forward-compatible: any in-flight `wait_for_node` / `wait_for_cache_task` retries queued in RabbitMQ at deploy time execute and return harmlessly (parent/cache will be settled by the time they run).
- `NodeTask` records for analyses in flight during deploy may carry no lease; the first beat sweep (or any user edit) re-dispatches them. Stale rows are cleaned by `delete_analysis_old_node_versions`.
- Add the `reschedule_stalled_analyses` beat entry when deploying; without it, lease-expiry reclamation will not fire (event-driven and entry-point triggers still work, but dead-worker self-healing depends on the periodic sweep).

---

## Follow-ups

- **#1553** — once this lands, `_node_launch_signature` / the per-layer lease round is the choke-point for fusing bulk node loads into CTE-shared SQL (create-from-template, force reload). Contingent on this issue; the seam does not exist in the current scheduler.
