import logging
import uuid
from collections import defaultdict
from datetime import timedelta
from typing import Optional

import celery
from auditlog.context import disable_auditlog
from celery.canvas import Signature, chain
from django.conf import settings
from django.db import transaction
from django.db.models import F, Q
from django.utils import timezone

from analysis.models import AnalysisEdge, AnalysisNode, NodeColors, NodeStatus, NodeTask
from analysis.models.nodes.analysis_node import LEASE_SECONDS, NodeCache, NodeVersion
from analysis.models.nodes.node_utils import get_nodes_by_id
from analysis.tasks.node_update_tasks import MAX_NODE_ATTEMPTS
from library.log_utils import log_traceback, report_message
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import JobsControl, ProcessingStatus


@celery.shared_task
def create_and_launch_analysis_tasks(analysis_id, run_async=True, max_nodes=None):
    """ Routed to scheduling_single_worker (the only caller of lease_ready_nodes).

        mocha dispatch_pending_actions(): lease the nodes that are ready right now, then fan
        out one task per leased node so analysis_workers pick them up in parallel. A node
        finishing re-triggers this dispatcher, which leases whatever just became ready.

        max_nodes bounds one dispatch (default ANALYSIS_NODE_DISPATCH_MAX_NODES_PER_ANALYSIS): a
        graph with hundreds of simultaneously-ready nodes goes out in waves, each wave's
        completions pulling through the next. Returns how many node tasks were launched. """
    if max_nodes is None:
        max_nodes = settings.ANALYSIS_NODE_DISPATCH_MAX_NODES_PER_ANALYSIS
    try:
        worker_id = f"dispatch:{create_and_launch_analysis_tasks.request.id or 'sync'}"
        leased = lease_ready_nodes(analysis_id, worker_id, max_nodes=max_nodes)
    except Exception:
        log_traceback()
        raise

    logging.info("create_and_launch_analysis_tasks(%s): leased %d node(s): %s",
                 analysis_id, len(leased), leased)
    launched = 0
    for node_id, version in leased:
        try:
            sig = _node_launch_signature(node_id, version)
        except AnalysisNode.DoesNotExist:
            # Version bumped (user edit) or node deleted between leasing and dispatch. The new
            # version is DIRTY and will be re-dispatched; skip so the other leased nodes still launch.
            continue
        launched += 1
        if run_async:
            sig.apply_async()
        else:
            result = sig.apply()
            if not result.successful():
                raise Exception(result.result)
    return launched


def _live_lease_q(now) -> Q:
    """ Node is being actively worked: its CURRENT version holds a lease that hasn't lapsed.
        Leases on superseded versions say nothing about the node's current state. """
    return Q(nodeversion__version=F("version"), nodeversion__nodetask__lease_expires__gte=now)


def _in_flight_node_count(now, limit) -> int:
    """ Nodes a worker is actively working (or about to), counted with a LIMIT - we only ever
        compare against the target, so stop as soon as we reach it. """
    qs = (AnalysisNode.objects
          .filter(status__in=NodeStatus.LOADING_STATUSES)
          .exclude(status=NodeStatus.DIRTY)
          .filter(_live_lease_q(now))
          .values_list("pk", flat=True)[:limit])
    return len(qs)


@celery.shared_task
def dispatch_analysis_backlog():
    """ Periodic backlog drain + dead-worker self-heal (issue #346), routed to
        scheduling_single_worker so it leases inline rather than fanning out a dispatch task per
        analysis.

        Analyses someone is working on dispatch themselves (update_analysis, and each node
        completion re-triggering its own analysis), so this only ever spends capacity that live
        work is not using: it tops the in-flight node count up to
        ANALYSIS_NODE_DISPATCH_BACKLOG_IN_FLIGHT_TARGET and exits immediately once that is met.
        Work waiting on a lapsed lease is not in flight, so reclaiming a dead worker's node still
        happens here.

        Candidates are deliberately over-inclusive (backoff windows and the reclaim / re-lease /
        terminal-fail decisions are all enforced in lease_ready_nodes) but bounded per tick. """
    if JobsControl.is_paused():
        return  # operational brake (e.g. crash safety auto-pause) - dispatch nothing
    now = timezone.now()
    target = settings.ANALYSIS_NODE_DISPATCH_BACKLOG_IN_FLIGHT_TARGET
    budget = target - _in_flight_node_count(now, target)
    if budget <= 0:
        logging.debug("dispatch_analysis_backlog: %d node(s) in flight - leaving the workers alone", target)
        return

    analysis_ids = (AnalysisNode.objects
                    .filter(status__in=NodeStatus.LOADING_STATUSES)
                    .exclude(_live_lease_q(now))
                    .order_by("analysis_id")
                    .values_list("analysis_id", flat=True)
                    .distinct()[:settings.ANALYSIS_NODE_DISPATCH_BACKLOG_MAX_ANALYSES])

    launched = 0
    for analysis_id in analysis_ids:
        launched += create_and_launch_analysis_tasks(analysis_id, max_nodes=budget - launched)
        if launched >= budget:
            break
    if launched:
        logging.info("dispatch_analysis_backlog: launched %d node task(s) of a %d budget", launched, budget)


def _load_graph(analysis_id):
    """ Returns (nodes_by_id, parents_by_child):
        - nodes_by_id: subclass instances keyed by pk (current statuses + cache-aware methods)
        - parents_by_child: {child_id: {parent_id, ...}} from AnalysisEdge """
    nodes_qs = AnalysisNode.objects.filter(analysis_id=analysis_id).select_subclasses()
    nodes_by_id = get_nodes_by_id(nodes_qs)
    parents_by_child = defaultdict(set)
    edges = AnalysisEdge.objects.filter(parent__analysis_id=analysis_id).values_list("parent_id", "child_id")
    for parent_id, child_id in edges:
        parents_by_child[child_id].add(parent_id)
    return nodes_by_id, parents_by_child


def _node_cache_target(node, force_cache=False) -> Optional[tuple]:
    """ The (node_id, version) of the NodeVersion whose NodeCache this node depends on, or
        None if the node uses no cache. Mirrors get_cache_task_args_set's target resolution
        WITHOUT creating anything. """
    if not (node.is_valid and (force_cache or node.use_cache)):
        return None
    if parent := node.get_unmodified_single_parent_node():
        return _node_cache_target(parent, force_cache=force_cache)
    return node.pk, node.version


def _node_cache_ready(node) -> bool:
    """ Hold a node back only while a shared NodeCache it depends on is actively being built by
        ANOTHER node (ProcessingStatus.PROCESSING); the builder's node_cache_task re-triggers the
        dispatcher on completion.

        A node that builds its own cache is never gated here - node_cache_task is chained ahead of
        its update task by _node_launch_signature, so gating it would deadlock a reclaimed builder.
        A cache that merely exists as CREATED (not yet, or no longer, being built) also does not
        block: the node falls back to a live query rather than hanging if the builder was lost. """
    target = _node_cache_target(node)
    if target is None or target == (node.pk, node.version):
        return True  # no cache dependency, or this node builds its own (chained, not gated)
    target_node_id, target_version = target
    node_cache = NodeCache.objects.filter(node_version__node_id=target_node_id,
                                          node_version__version=target_version) \
        .select_related("variant_collection").first()
    if node_cache is None:
        return True  # not created yet - this node (or a same-round sibling) will build it
    return node_cache.variant_collection.status != ProcessingStatus.PROCESSING


def _node_ready_to_lease(node, parents_by_child, nodes_by_id) -> bool:
    """ mocha dependencies_satisfied() equivalent.

        All parents must be settled (in READY_STATUSES - READY or any error), i.e. none still
        loading. ERROR parents count as 'ready' so the child runs and fails fast with
        ERROR_WITH_PARENT (matches issue comment). A required NodeCache being built by another
        node must also have finished. """
    for parent_id in parents_by_child.get(node.id, ()):
        parent = nodes_by_id.get(parent_id)
        if parent is None or NodeStatus.is_loading(parent.status):
            return False
    return _node_cache_ready(node)


def lease_ready_nodes(analysis_id, worker_id, lease_seconds=LEASE_SECONDS, max_nodes=None):
    """ mocha lease_action() + dispatch, batched per analysis. SINGLE-WORKER ONLY (called only
        from create_and_launch_analysis_tasks on scheduling_single_worker).

        Claims fresh work, reclaims expired leases, and gives up (terminal fail) past
        MAX_NODE_ATTEMPTS. max_nodes caps how many are claimed - the rest stay DIRTY for the next
        dispatch. Returns the (node_id, version) tuples now owned by this run. """
    if max_nodes is not None and max_nodes <= 0:
        return []
    if JobsControl.is_paused():
        return []  # operational brake (e.g. crash safety auto-pause) - lease nothing
    now = timezone.now()
    leased = []
    leased_cache_targets = set()
    with transaction.atomic():
        # Coarse base-table filter (no select_subclasses - keeps select_for_update simple).
        # Candidates are nodes that are not yet settled (DIRTY/QUEUED/LOADING_CACHE/LOADING).
        # skip_locked => a racing dispatcher simply skips rows we hold.
        # select_related(None) clears the manager's nullable LEFT JOINs (annotation_version etc.)
        # which Postgres refuses to lock ("FOR UPDATE cannot be applied to the nullable side...").
        candidates = list(
            AnalysisNode.objects
            .select_related(None)
            .select_for_update(skip_locked=True)
            .filter(analysis_id=analysis_id, status__in=NodeStatus.LOADING_STATUSES)
            .order_by("pk")
        )
        nodes_by_id, parents_by_child = _load_graph(analysis_id)
        tasks_by_node = {}  # current-version lease rows only
        for nt in NodeTask.objects.filter(node_version__node__analysis_id=analysis_id).select_related("node_version"):
            node = nodes_by_id.get(nt.node_version.node_id)
            if node and node.version == nt.node_version.version:
                tasks_by_node[nt.node_version.node_id] = nt

        for candidate in candidates:
            node = nodes_by_id.get(candidate.pk)
            if node is None:
                continue
            node_task = tasks_by_node.get(node.pk)
            lease_live = bool(node_task and node_task.lease_expires and node_task.lease_expires >= now)
            if node.status != NodeStatus.DIRTY and lease_live:
                continue  # someone is actively working it, lease still valid -> leave alone

            # Two kinds of abandoned (non-DIRTY, lease lapsed) work, told apart by the status the
            # worker durably committed before it vanished:
            #  - LOADING / LOADING_CACHE: a worker started executing THIS node and never returned
            #    (hard machine lockup, OOM-killed process, SIGKILL). The node is the prime suspect
            #    for having killed the worker, so perma-fail it rather than re-dispatching it and
            #    crashing the next worker too. (An in-process OOM is caught in update_node_task and
            #    already fails cleanly without reaching here; this is the backstop for hard death.)
            #  - QUEUED (falls through below): the worker died before it began this node (deploy /
            #    restart / spot reclaim) - the node never ran, so re-leasing it is correct self-heal.
            if node.status in (NodeStatus.LOADING, NodeStatus.LOADING_CACHE):
                msg = (f"Worker lost while loading node {node.pk}/{node.version} "
                       f"(analysis {analysis_id}) - possible out-of-memory or machine lockup. "
                       f"Permanently failed to avoid repeatedly crashing workers; restart manually if needed.")
                _fail_node(node, msg)
                report_message(msg, level='error')
                continue

            attempts = node_task.attempt_count if node_task else 0
            if attempts >= MAX_NODE_ATTEMPTS:
                # Burned all attempts (transient errors and/or dead-worker reclaims) -> give up.
                _fail_node(node, "Node task did not complete (worker lost); marked failed.")
                continue
            if node_task and node_task.run_after and node_task.run_after > now:
                continue  # backing off - not leasable yet
            if not _node_ready_to_lease(node, parents_by_child, nodes_by_id):
                continue  # parents still loading / required cache still building

            cache_target = _node_cache_target(node)
            is_own_cache = cache_target == (node.pk, node.version)
            if cache_target is not None and not is_own_cache and cache_target in leased_cache_targets:
                continue  # a node leased this round is building the shared cache - wait for it

            node_version, _ = NodeVersion.objects.get_or_create(node_id=node.pk, version=node.version)
            NodeTask.objects.update_or_create(
                node_version=node_version,
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
            if cache_target is not None:
                leased_cache_targets.add(cache_target)
            if max_nodes is not None and len(leased) >= max_nodes:
                break
    return leased


def _fail_node(node, message):
    with disable_auditlog():
        AnalysisNode.objects.filter(pk=node.pk, version=node.version).update(
            status=NodeStatus.ERROR, shadow_color=NodeColors.ERROR, errors=message)
    NodeTask.objects.filter(node_version__node_id=node.pk, node_version__version=node.version) \
        .update(lease_expires=None, leased_by=None)


def _node_launch_signature(node_id, version):
    """ Per-node replacement for the old group/chord building. For an ordinary node returns its
        update task; for a node that must build its OWN cache first returns the two-task chain
        [node_cache_task, update_node_task] (a single-node chain - no cross-node dependency, so
        no deadlock risk). Shared caches are deduped by NodeCache.get_or_create_for_node + the
        in-flight-cache dependency gate, not by a wait task. """
    node = AnalysisNode.objects.get_subclass(pk=node_id, version=version)
    update_job = node.get_update_task()
    cache_jobs = [Signature(t, args=a, immutable=True)
                  for (t, a) in node.get_cache_task_args_set() if t]
    if cache_jobs:
        return chain(*cache_jobs, update_job)  # same node: cache then update
    return update_job


@celery.shared_task
def populate_clingen_alleles_from_analysis_node(node_id, max_variants=0):
    node = AnalysisNode.objects.get_subclass(pk=node_id)
    variants_qs = node.get_queryset()
    if max_variants:
        # If we're going to limit it, make sure we only get useful ones....
        variants_qs = variants_qs.filter(variantallele__allele__clingen_allele__isnull=True)[:max_variants]
    populate_clingen_alleles_for_variants(node.analysis.genome_build, variants_qs)
