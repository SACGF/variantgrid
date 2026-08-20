# Analysis Node Scheduling Plan — Issue #346

## Problem

When nodes are added to an analysis while other nodes are already running, the scheduler inserts `wait_for_node` tasks into Celery chains to hold children back until parents finish. The original implementation used `sleep()` inside these tasks, **occupying a Celery worker** the entire wait duration. If enough accumulated, the parent nodes could never get a worker — a deadlock.

The VG.com crash log showed exactly this:

```
wait_for_node   Waiting on parent node 78301 which is QUEUED - waiting for 5 secs...
wait_for_node   Waiting on parent node 78300 which is QUEUED - waiting for 5 secs...
wait_for_node   Waiting on parent node 78300 which is QUEUED - waiting for 5 secs...
```

**Root cause**: The scheduler builds a static chain that embeds blocking wait-tasks to represent a dynamic dependency — fundamentally at odds with a DAG that users can modify mid-run.

---

## Stop-gap Already Applied (variantgrid_com#90)

`wait_for_node` and `wait_for_cache_task` were converted from `sleep()`-based polling to `self.retry(countdown=N)` (commit in variantgrid_com#90). The tasks now re-queue themselves with a delay and **release the worker immediately** between checks.

```python
# wait_for_node (after stop-gap)
@celery.shared_task(bind=True)
def wait_for_node(self, node_id):
    ...
    raise self.retry(countdown=sleep_time, max_retries=len(TIME_BETWEEN_CHECKS))

# wait_for_cache_task (after stop-gap)
@celery.shared_task(bind=True)
def wait_for_cache_task(self, node_cache_id):
    ...
    raise self.retry(countdown=1, max_retries=MAX_CHECKS)
```

This removes the crash vector (worker starvation from blocking sleeps). However, there is one remaining blocking call in `wait_for_node`: when a `NodeTask.celery_task` exists, it still calls `wait_for_task()` which does `AsyncResult.get(timeout=300s)` — blocking the worker for up to 5 minutes:

```python
# Still blocks a worker for up to 5 minutes:
node_task = NodeTask.objects.get(node=node, version=node.version)
if node_task.celery_task:
    wait_for_task(node_task.celery_task)  # ← result.get(timeout=MINUTE_SECS * 5)
    return
```

This path is taken when the parent node has already started executing (has a real Celery task ID). It is a lower-risk scenario than the QUEUED case (fewer workers blocked, shorter expected wait) but is still a worker-occupying pattern.

**The architectural issue remains**: `wait_for_node` tasks still occupy queue slots and execute repeatedly while waiting. The proper fix is to never schedule tasks with unmet dependencies.

---

## Solution: Post-Completion Re-Scheduling

**Core principle**: Never dispatch a task whose dependencies are unmet. Instead, each time a node finishes (success or error), trigger the scheduler to pick up any newly-unblocked work.

This replaces the blocking `wait_for_node` / `wait_for_cache_task` pattern with an event-driven loop:

```
node completes → trigger create_and_launch_analysis_tasks(analysis_id)
                → scheduler finds nodes whose parents just became READY
                → schedules them
                → they complete → trigger scheduler again
                → ...
```

**Why this is safe** — `create_and_launch_analysis_tasks` already has a distributed lock:
```python
# NodeTask.unique_together = ("node", "version")
NodeTask.objects.bulk_create(node_task_records, ignore_conflicts=True)
node_tasks = NodeTask.objects.filter(analysis_update_uuid=analysis_update_uuid)
node_versions_to_update = dict(node_tasks.values_list("node_id", "version"))
```
Two scheduler runs racing for the same node: only one inserts successfully; the other gets an empty `node_versions_to_update` and returns quickly. Combined with `scheduling_single_worker` (single concurrency), multiple triggers for the same analysis serialize with minimal overhead.

---

## Files Changed

| File | Change |
|------|--------|
| `analysis/tasks/analysis_update_tasks.py` | Build `parent_value_data` earlier; filter node locks to schedulable-only; remove `wait_for_node` from chains; remove `_get_node_dependencies` helper |
| `analysis/tasks/node_update_tasks.py` | `update_node_task` and `node_cache_task` trigger re-scheduling on completion |

---

## Detailed Changes

### 1. `_can_schedule_node` — new helper

A node can run now only if all its parents are in `READY_STATUSES` (not loading):

```python
# analysis/tasks/analysis_update_tasks.py

def _can_schedule_node(node_id, nodes_by_id, parent_value_data):
    """Returns True if all parents are settled (READY_STATUSES) — not loading."""
    for parent_id in parent_value_data.get(node_id, set()):
        parent = nodes_by_id.get(parent_id)
        if parent and NodeStatus.is_loading(parent.status):
            return False
    return True
```

`NodeStatus.LOADING_STATUSES = [DIRTY, QUEUED, LOADING_CACHE, LOADING]` — all four prevent scheduling children.

---

### 2. `_get_analysis_update_tasks` — build `parent_value_data` earlier, filter locks

The `parent_value_data` dict (edges in the sub-graph) is currently built *after* the NodeTask lock loop. Move it before so `_can_schedule_node` can use it:

```python
# analysis/tasks/analysis_update_tasks.py  _get_analysis_update_tasks()

for connected_components in nx.weakly_connected_components(all_nodes_graph):
    sub_graph = all_nodes_graph.subgraph(connected_components)
    sub_graph_node_ids = list(sub_graph)

    # Build parent_value_data FIRST (moved up from below)
    parent_value_data = defaultdict(set)
    for parent, child_list in nx.to_dict_of_lists(sub_graph).items():
        for child_node_id in child_list:
            parent_value_data[child_node_id].add(parent)

    analysis_update_uuid = uuid.uuid4()
    node_task_records = []
    for node_id in sub_graph_node_ids:
        node = nodes_by_id[node_id]
        # OLD: if node.status == NodeStatus.DIRTY:
        # NEW: only lock nodes whose parents are already settled
        if node.status == NodeStatus.DIRTY and _can_schedule_node(node_id, nodes_by_id, parent_value_data):
            node_task = NodeTask(node_id=node_id, version=node.version, analysis_update_uuid=analysis_update_uuid)
            node_task_records.append(node_task)

    # ... rest of locking / toposort unchanged ...

    groups = []
    if node_versions_to_update:
        # parent_value_data already built above — no need to rebuild
        topo_sorted = get_toposorted_nodes_from_parent_value_data(nodes_by_id, parent_value_data)

        all_cache_jobs = set()
        for grp in topo_sorted:
            group_cache_jobs = _add_jobs_for_group(node_versions_to_update, grp, groups, all_cache_jobs)
            all_cache_jobs.update(group_cache_jobs)
        # ... QUEUED update unchanged ...
```

**Key difference**: dirty nodes with loading parents are not locked into `NodeTask`. When the parent completes and the scheduler re-runs, those children will be seen as DIRTY with no existing `NodeTask`, so they can be locked and scheduled then.

---

### 3. `_add_jobs_for_group` — remove `wait_for_node` branch

**Old signature**: `_add_jobs_for_group(nodes_to_update, dependencies, grp, groups, existing_cache_jobs)`
**New signature**: `_add_jobs_for_group(nodes_to_update, grp, groups, existing_cache_jobs)`
(drop the `dependencies` parameter — it only existed to support `wait_for_node`)

```python
# analysis/tasks/analysis_update_tasks.py

def _add_jobs_for_group(nodes_to_update, grp, groups, existing_cache_jobs) -> set:
    cache_jobs = set()
    after_jobs = []
    jobs = []

    for node in grp:
        if node.pk in nodes_to_update:
            update_job = node.get_update_task()
            task_args_objs_set = node.get_cache_task_args_set()
            new_cache_jobs = task_args_objs_set - existing_cache_jobs

            # NEW: if cache is already being built elsewhere, skip this node.
            # node_cache_task will trigger the scheduler when done.
            wait_jobs = {(t, a) for t, a in new_cache_jobs if t == AnalysisNode.WAIT_FOR_CACHE_TASK}
            if wait_jobs:
                # Don't schedule: release the lock so next scheduler run can re-acquire it
                nodes_to_update.pop(node.pk, None)
                continue

            build_jobs = new_cache_jobs - wait_jobs  # NODE_CACHE_TASK entries
            if build_jobs:
                cache_jobs.update(build_jobs)
                after_jobs.append(update_job)
            else:
                jobs.append(update_job)

        # REMOVED: elif node.pk in dependencies / wait_for_node branch

    if cache_jobs:
        for task, args in cache_jobs:
            if task:
                jobs.append(Signature(task, args=args, immutable=True))

    groups.append(jobs)
    groups.append(after_jobs)
    return cache_jobs
```

> **Note on `WAIT_FOR_CACHE_TASK` case**: when `get_cache_task_args_set()` returns `WAIT_FOR_CACHE_TASK`, the Venn/intersection node's cache is already being built by another task run. Rather than inserting a blocking wait task, we drop the lock on this child node (`nodes_to_update.pop`) so the upcoming `node_cache_task` completion can re-schedule it. The child stays DIRTY; the running cache task's post-completion trigger will pick it up.

Also remove the now-unused helpers:

```python
# DELETE these functions:
def _get_ancestor_set(node_id, parent_value_data): ...
def _get_node_dependencies(nodes_by_id, parent_value_data): ...
```

And the import at the top:

```python
# DELETE:
from analysis.tasks.node_update_tasks import wait_for_node
```

---

### 4. `update_node_task` — trigger re-scheduling on completion

After every outcome (success, error, cancellation), fire the scheduler for the analysis. Also fire a delayed copy as a safety net for the race where the parent finishes just as the scheduler reads node states.

```python
# analysis/tasks/node_update_tasks.py

@celery.shared_task(base=AbortableTask)
def update_node_task(node_id, version):
    with disable_auditlog():
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
        except AnalysisNode.DoesNotExist:
            return  # stale / deleted

        analysis_id = node.analysis_id  # capture before load() or deletion

        errors = None
        node_errors = node.get_errors()
        if not node_errors:
            try:
                node.set_node_task_and_status(update_node_task.request.id, NodeStatus.LOADING)
                node.load()
                if node.shadow_color == NodeColors.ERROR and node.is_valid:
                    node.update(shadow_color=None)
                return  # success path
            except NodeOutOfDateException:
                logging.warning("Node %d/%d out of date - exiting", node.pk, node.version)
                return  # version bumped — scheduler will re-run for the new version
            except OperationalError:
                status = NodeStatus.CANCELLED
            except NodeConfigurationException:
                status = NodeStatus.ERROR_CONFIGURATION
            except NodeParentErrorsException:
                status = NodeStatus.ERROR_WITH_PARENT
            except Exception:
                errors = get_traceback()
                status = NodeStatus.ERROR
        else:
            status = AnalysisNode.get_status_from_errors(node_errors)

        try:
            shadow_color = NodeColors.ERROR if NodeStatus.is_error(status) else None
            node.update(status=status, errors=errors, shadow_color=shadow_color)
        except (IntegrityError, NodeOutOfDateException):
            pass

    # NEW: trigger scheduler so children blocked on this node can now run.
    # The immediate trigger handles the common path; the delayed one is a safety
    # net for the race where this task completes just before the scheduler reads
    # node statuses (comment from issue #346).
    _trigger_rescheduling(analysis_id)


def _trigger_rescheduling(analysis_id):
    """Schedule two re-scheduling attempts: one immediate and one with a short delay.
    Both are no-ops if there is nothing to do (NodeTask lock prevents double-dispatch).
    Runs in scheduling_single_worker so multiple triggers for the same analysis serialize."""
    from analysis.tasks.analysis_update_tasks import create_and_launch_analysis_tasks
    create_and_launch_analysis_tasks.apply_async(args=[analysis_id])
    create_and_launch_analysis_tasks.apply_async(args=[analysis_id], countdown=3)
```

> The `return` on the success path means `_trigger_rescheduling` is not called. Move it outside the `try/except` block, or restructure with a `finally`. The simplest restructure: call `_trigger_rescheduling(analysis_id)` **after** the function body in a `finally` clause (but only if `analysis_id` was resolved).

Revised structure using a flag:

```python
@celery.shared_task(base=AbortableTask)
def update_node_task(node_id, version):
    analysis_id = None
    with disable_auditlog():
        try:
            node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
            analysis_id = node.analysis_id
        except AnalysisNode.DoesNotExist:
            return

        # ... (all existing logic, including the early returns for NodeOutOfDateException) ...

    if analysis_id is not None:
        _trigger_rescheduling(analysis_id)
```

For `NodeOutOfDateException` early returns: the version was bumped (user edited the node), which triggers a new `reload_analysis_nodes` → `update_analysis` → fresh scheduling pass. No need to trigger manually.

---

### 5. `node_cache_task` — trigger re-scheduling on completion

```python
# analysis/tasks/node_update_tasks.py

@celery.shared_task
def node_cache_task(node_id, version):
    node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
    analysis_id = node.analysis_id  # capture before any exception
    node_cache = NodeCache.objects.get(node_version=node.node_version)
    variant_collection = node_cache.variant_collection

    if variant_collection.status != ProcessingStatus.CREATED:
        _trigger_rescheduling(analysis_id)
        return

    if not (node.is_valid and node.modifies_parents()):
        variant_collection.status = ProcessingStatus.SKIPPED
        variant_collection.save()
        _trigger_rescheduling(analysis_id)
        return

    # ... (rest of existing logic unchanged) ...

    variant_collection.status = processing_status
    variant_collection.save()

    if node.status != NodeStatus.READY:
        with disable_auditlog():
            node.status = status
            node.save()

    # NEW: unblock children that were waiting for this cache
    _trigger_rescheduling(analysis_id)
```

After `node_cache_task` completes (either SKIPPED, SUCCESS, or ERROR), the scheduler re-runs. If the corresponding `update_node_task` ran after (in the chain), the node will be READY at that point and children whose only blocker was the cache can now be scheduled.

---

## Handling the `wait_for_cache_task` node lock release

The `WAIT_FOR_CACHE_TASK` skip in `_add_jobs_for_group` removes the node from `nodes_to_update`. But `NodeTask` records were already bulk-inserted before `_add_jobs_for_group` runs. We need to clean up the lock for skipped nodes so the next scheduler run can re-acquire it.

Two options:

**Option A** — Delete the NodeTask for skipped nodes:
```python
if wait_jobs:
    NodeTask.objects.filter(node_id=node.pk, version=node.version,
                            analysis_update_uuid=...).delete()
    nodes_to_update.pop(node.pk, None)
    continue
```
This requires passing `analysis_update_uuid` into `_add_jobs_for_group`.

**Option B** — Don't lock nodes with pending caches in the first place.

Option B is cleaner. Add a cache-readiness check to the lock filter in `_get_analysis_update_tasks`:

```python
def _can_schedule_node(node_id, nodes_by_id, parent_value_data):
    """Returns True if all parents are settled AND no required caches are still processing."""
    for parent_id in parent_value_data.get(node_id, set()):
        parent = nodes_by_id.get(parent_id)
        if parent and NodeStatus.is_loading(parent.status):
            return False
    return True  # cache readiness checked lazily in _add_jobs_for_group if needed
```

For now, implementing Option A is simpler (minimal change to `_can_schedule_node`). The `WAIT_FOR_CACHE_TASK` case is rare (requires a second scheduler run while the cache task is mid-flight) and the lock deletion just unblocks the next pass.

---

## Race Condition Analysis

**Scenario**: Parent node A finishes just *before* the scheduler runs (so the scheduler's `nodes_by_id` snapshot shows A as LOADING, not READY).

Old code: child B had `wait_for_node(A)` in its chain — it would still work, just poll.

New code: scheduler sees A as loading, skips B. Then the scheduler returns. Nobody re-triggers.

**Mitigation**: `_trigger_rescheduling` always dispatches **two** triggers:
```python
create_and_launch_analysis_tasks.apply_async(args=[analysis_id])           # immediate
create_and_launch_analysis_tasks.apply_async(args=[analysis_id], countdown=3)  # 3s delay
```

The sequence:
1. A starts loading (status = LOADING)
2. Scheduler 1 runs: A is loading → B not scheduled
3. A completes: status = READY; triggers Scheduler 2 (immediate) + Scheduler 3 (3s delay)
4. Scheduler 2 runs: A is READY → B scheduled ✓

If Scheduler 2 runs while A's DB write is not yet committed (unlikely with `scheduling_single_worker` serialisation but possible in theory), Scheduler 3 fires 3 seconds later and catches it.

The delayed trigger is the "safety net" mentioned in the issue comments:
> *run another job with a small delay to guarantee something executes and performs the scheduling. Since it runs through the single worker thread we can exit here if already done*

---

## State Transitions (new model)

```
User edits analysis
    └─> reload_analysis_nodes()
            └─> all nodes → DIRTY
            └─> update_analysis(analysis_id)
                    └─> create_and_launch_analysis_tasks (scheduling_single_worker)
                            └─> for each DIRTY node with all-READY parents:
                                    lock NodeTask → QUEUED
                                    dispatch update_node_task chain
                                                │
                                                ▼
                                    update_node_task → LOADING → READY/ERROR
                                                │
                                                ▼
                                    _trigger_rescheduling()
                                                │
                                                ▼
                                    create_and_launch_analysis_tasks (scheduling_single_worker)
                                            └─> find newly-unblocked DIRTY nodes → schedule them
                                            └─> (if nothing to do → fast exit)
```

---

## What is NOT Changed

- `wait_for_node` task: kept for `analysis_grid_export_tasks._wait_for_output_node()`, which calls it synchronously as a blocking export gate (not in a chain). Also kept as the implementation but it will no longer be added to scheduling chains.
- The stop-gap `self.retry()` conversion in `wait_for_node` / `wait_for_cache_task` is preserved — those tasks still exist and work correctly for any edge cases; they just won't be inserted into new chains.
- `wait_for_task()` / `AsyncResult.get()` in `wait_for_node`: this blocking call is a lower-priority concern once the task is no longer added to chains. It remains for the synchronous export use case.
- `NodeTask.unique_together` lock: unchanged, still the distributed-lock mechanism.
- `node_cache_task` chaining *within* a node's own update (e.g. `[cache_task, update_task]` for Venn/Intersection nodes): the cache task still runs in the chain before the update task, because both are for the same node. Children of the Venn node are not included in this chain — they're scheduled on the next pass after Venn's `update_node_task` fires `_trigger_rescheduling`.
- `reload_analysis_nodes` / `update_analysis`: unchanged entry points.

---

## Testing

1. **Unit test — `_can_schedule_node`**: dirty node with LOADING parent → False; dirty node with READY parent → True; dirty node with ERROR parent → True (error is a READY_STATUS, child should run and fail with `ERROR_WITH_PARENT`).

2. **Integration test — cascade scheduling**: create a 3-node chain (A→B→C). Mark A as DIRTY. Run `create_and_launch_analysis_tasks`. Verify only A is dispatched. Simulate A completing (set A to READY). Run scheduler again. Verify B is dispatched, C is not. Simulate B completing. Verify C is dispatched.

3. **Integration test — concurrent schedule calls**: two calls to `create_and_launch_analysis_tasks` for the same analysis at the same time. Verify each node is dispatched exactly once (NodeTask lock).

4. **Regression test — no `wait_for_node` in chains**: after any call to `_get_analysis_update_tasks`, verify no dispatched task is of type `wait_for_node`.

5. **Stress test (manual)**: load an auto-analysis scenario (as described in the issue). Previously caused crash. Verify all nodes reach READY or ERROR without worker starvation.

---

## Migration Notes

- No DB migrations needed.
- The stop-gap (`self.retry()` conversion) is already deployed to variantgrid_com. This full fix is safe to deploy over the top: in-flight `wait_for_node` / `wait_for_cache_task` retries already queued in RabbitMQ will execute and return harmlessly (parent will already be READY by the time they run, causing an immediate early return).
- `NodeTask` records for in-flight analyses during deployment may be orphaned. Safe to let them time out naturally; stale records are cleaned up by `delete_analysis_old_node_versions`.
