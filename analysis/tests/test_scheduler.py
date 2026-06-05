"""
Tests for analysis node scheduling logic.

Covers NodeStatus categorisation and the state-driven leased scheduler
(lease_ready_nodes + the dependency gate) introduced by issue #346.
"""
import inspect
import uuid
from datetime import timedelta
from unittest import mock

from celery.canvas import Signature, _chain
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.models import AllVariantsNode, NodeTask
from analysis.models.enums import NodeStatus, SetOperations
from analysis.models.nodes.analysis_node import NodeCache
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.models.nodes.filters.venn_node import VennNode, VennNodeCache, venn_cache_count
from analysis.tasks import analysis_update_tasks
from analysis.tasks.analysis_update_tasks import lease_ready_nodes, _node_launch_signature, \
    _node_ready_to_lease, create_and_launch_analysis_tasks
from analysis.tasks.node_update_tasks import next_backoff, _backoff_node, reschedule_stalled_analyses, \
    update_node_task, node_cache_task, MAX_NODE_ATTEMPTS
from analysis.tests.utils import AnalysisSetupMixin
from snpdb.models import ProcessingStatus, VariantCollection


class TestNodeStatusCoverage(TestCase):
    """Every NodeStatus value must be in exactly one of LOADING_STATUSES or READY_STATUSES."""

    def test_no_status_in_both_buckets(self):
        loading = set(NodeStatus.LOADING_STATUSES)
        ready = set(NodeStatus.READY_STATUSES)
        overlap = loading & ready
        self.assertEqual(overlap, set(), f"Status(es) appear in both LOADING and READY: {overlap}")

    def test_all_statuses_covered(self):
        all_values = {choice[0] for choice in NodeStatus.choices}
        loading = set(NodeStatus.LOADING_STATUSES)
        ready = set(NodeStatus.READY_STATUSES)
        uncategorised = all_values - loading - ready
        self.assertEqual(uncategorised, set(), f"Status(es) in neither bucket: {uncategorised}")

    def test_dirty_is_loading(self):
        self.assertTrue(NodeStatus.is_loading(NodeStatus.DIRTY))

    def test_queued_is_loading(self):
        self.assertTrue(NodeStatus.is_loading(NodeStatus.QUEUED))

    def test_loading_cache_is_loading(self):
        self.assertTrue(NodeStatus.is_loading(NodeStatus.LOADING_CACHE))

    def test_loading_is_loading(self):
        self.assertTrue(NodeStatus.is_loading(NodeStatus.LOADING))

    def test_ready_is_ready(self):
        self.assertTrue(NodeStatus.is_ready(NodeStatus.READY))

    def test_error_is_ready(self):
        # Error states are "done" — not loading — so children can proceed (and fail with ERROR_WITH_PARENT)
        self.assertTrue(NodeStatus.is_ready(NodeStatus.ERROR))
        self.assertTrue(NodeStatus.is_ready(NodeStatus.ERROR_CONFIGURATION))
        self.assertTrue(NodeStatus.is_ready(NodeStatus.ERROR_WITH_PARENT))
        self.assertTrue(NodeStatus.is_ready(NodeStatus.CANCELLED))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestSchedulerTasks(AnalysisSetupMixin, TestCase):
    """Integration tests for lease_ready_nodes() — the single place that claims work.

    Each test creates nodes fresh (they are rolled back after each test by
    Django's TestCase transaction handling). Statuses are forced via .update()
    to bypass AnalysisNode.save() version-bump logic.
    """

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _force_status(self, node, status):
        """Set node status directly in DB, bypassing save() version-bump logic."""
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def _lease(self):
        """Lease ready nodes for the analysis (mirrors create_and_launch's claim step)."""
        return lease_ready_nodes(self.analysis.pk, "test-worker")

    def _make_chain(self, parent_status, child_status=NodeStatus.DIRTY):
        """Return (parent AllVariantsNode, child GeneListNode) with an edge between them."""
        parent = AllVariantsNode.objects.create(analysis=self.analysis)
        child = GeneListNode.objects.create(analysis=self.analysis)
        child.add_parent(parent)
        self._force_status(parent, parent_status)
        self._force_status(child, child_status)
        return parent, child

    # ------------------------------------------------------------------
    # Basic scheduling invariants
    # ------------------------------------------------------------------

    def test_no_nodes_leases_nothing(self):
        self.assertEqual(self._lease(), [])

    def test_no_dirty_nodes_leases_nothing(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        self._force_status(node, NodeStatus.READY)
        self.assertEqual(self._lease(), [])

    def test_single_dirty_node_is_leased(self):
        AllVariantsNode.objects.create(analysis=self.analysis)  # default status=DIRTY
        self.assertEqual(len(self._lease()), 1)

    def test_single_dirty_node_becomes_queued(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        self._lease()
        node.refresh_from_db()
        self.assertEqual(node.status, NodeStatus.QUEUED)

    def test_single_dirty_node_gets_node_task_lock(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        self._lease()
        self.assertTrue(NodeTask.objects.filter(node_version__node=node).exists())

    # ------------------------------------------------------------------
    # Double-dispatch: the live lease prevents scheduling the same node twice
    # ------------------------------------------------------------------

    def test_double_dispatch_creates_only_one_node_task(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        self._lease()
        self._lease()
        self.assertEqual(NodeTask.objects.filter(node_version__node=node).count(), 1)

    def test_double_dispatch_second_call_leases_nothing(self):
        AllVariantsNode.objects.create(analysis=self.analysis)
        self._lease()
        self.assertEqual(self._lease(), [])

    # ------------------------------------------------------------------
    # DAG: parent-child scheduling
    # ------------------------------------------------------------------

    def test_ready_parent_allows_child_to_be_scheduled(self):
        """Child whose only parent is READY should be scheduled."""
        _, child = self._make_chain(parent_status=NodeStatus.READY)
        self._lease()
        self.assertTrue(NodeTask.objects.filter(node_version__node=child).exists())

    def test_error_parent_allows_child_to_be_scheduled(self):
        """Error states are READY_STATUSES; child should still be scheduled
        (it will fail with ERROR_WITH_PARENT when it actually runs)."""
        _, child = self._make_chain(parent_status=NodeStatus.ERROR)
        self._lease()
        self.assertTrue(NodeTask.objects.filter(node_version__node=child).exists())

    def test_loading_parent_child_not_locked(self):
        """Child should NOT be scheduled (locked via NodeTask) when its parent
        is still loading.  The scheduler skips it and the post-completion trigger
        picks it up once the parent reaches READY (issue #346 fix)."""
        _, child = self._make_chain(parent_status=NodeStatus.LOADING)
        self._lease()
        self.assertFalse(
            NodeTask.objects.filter(node_version__node=child).exists(),
            "Child was locked (NodeTask created) despite parent being LOADING — "
            "wait_for_node pattern is still active",
        )

    def test_queued_parent_child_not_locked(self):
        """QUEUED is a LOADING_STATUS: child must not be scheduled.
        This is the exact scenario from the VG.com crash log."""
        _, child = self._make_chain(parent_status=NodeStatus.QUEUED)
        self._lease()
        self.assertFalse(
            NodeTask.objects.filter(node_version__node=child).exists(),
            "Child was locked despite parent being QUEUED",
        )


WAIT_FOR_NODE_TASK = "analysis.tasks.node_update_tasks.wait_for_node"


def _signature_task_names(sig):
    """Flatten a launch signature (Signature or chain) into its task names."""
    if isinstance(sig, _chain):
        return [t.task for t in sig.tasks]
    return [sig.task]


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestNodeReadyToLease(AnalysisSetupMixin, TestCase):
    """Unit tests for the dependency gate _node_ready_to_lease() (mocha dependencies_satisfied)."""

    def _force_status(self, node, status):
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def _graph(self):
        return analysis_update_tasks._load_graph(self.analysis.pk)

    def _make_chain(self, parent_status):
        parent = AllVariantsNode.objects.create(analysis=self.analysis)
        child = GeneListNode.objects.create(analysis=self.analysis)
        child.add_parent(parent)
        self._force_status(parent, parent_status)
        return parent, child

    def test_parent_loading_blocks_child(self):
        _, child = self._make_chain(NodeStatus.LOADING)
        nodes_by_id, parents_by_child = self._graph()
        self.assertFalse(_node_ready_to_lease(nodes_by_id[child.pk], parents_by_child, nodes_by_id))

    def test_parent_queued_blocks_child(self):
        _, child = self._make_chain(NodeStatus.QUEUED)
        nodes_by_id, parents_by_child = self._graph()
        self.assertFalse(_node_ready_to_lease(nodes_by_id[child.pk], parents_by_child, nodes_by_id))

    def test_parent_ready_allows_child(self):
        _, child = self._make_chain(NodeStatus.READY)
        nodes_by_id, parents_by_child = self._graph()
        self.assertTrue(_node_ready_to_lease(nodes_by_id[child.pk], parents_by_child, nodes_by_id))

    def test_parent_error_allows_child(self):
        """ERROR parents are settled — child runs and fails fast with ERROR_WITH_PARENT."""
        _, child = self._make_chain(NodeStatus.ERROR)
        nodes_by_id, parents_by_child = self._graph()
        self.assertTrue(_node_ready_to_lease(nodes_by_id[child.pk], parents_by_child, nodes_by_id))

    # --- cache dependency gate (replaces WAIT_FOR_CACHE_TASK) ---

    def _make_cache_owner(self, status):
        """A node owning a NodeCache whose variant_collection has the given ProcessingStatus."""
        owner = AllVariantsNode.objects.create(analysis=self.analysis)
        node_cache, _ = NodeCache.get_or_create_for_node(owner)
        node_cache.variant_collection.status = status
        node_cache.variant_collection.save()
        return owner

    def test_in_flight_shared_cache_blocks_consumer(self):
        """A consumer waits while ANOTHER node's shared cache is actively building (PROCESSING)."""
        owner = self._make_cache_owner(ProcessingStatus.PROCESSING)
        consumer = AllVariantsNode.objects.create(analysis=self.analysis)
        nodes_by_id, parents_by_child = self._graph()
        with mock.patch.object(analysis_update_tasks, "_node_cache_target",
                               return_value=(owner.pk, owner.version)):
            self.assertFalse(_node_ready_to_lease(nodes_by_id[consumer.pk], parents_by_child, nodes_by_id))

    def test_finished_shared_cache_allows_consumer(self):
        for status in (ProcessingStatus.SUCCESS, ProcessingStatus.SKIPPED):
            owner = self._make_cache_owner(status)
            consumer = AllVariantsNode.objects.create(analysis=self.analysis)
            nodes_by_id, parents_by_child = self._graph()
            with mock.patch.object(analysis_update_tasks, "_node_cache_target",
                                   return_value=(owner.pk, owner.version)):
                self.assertTrue(_node_ready_to_lease(nodes_by_id[consumer.pk], parents_by_child, nodes_by_id),
                                f"shared cache status {status} should not block")

    def test_node_building_own_cache_is_never_gated(self):
        """Liveness: a node that builds its own cache must not be blocked by that cache's status —
        otherwise a reclaimed builder (worker died mid-build) would hang forever."""
        owner = self._make_cache_owner(ProcessingStatus.PROCESSING)
        nodes_by_id, parents_by_child = self._graph()
        with mock.patch.object(analysis_update_tasks, "_node_cache_target",
                               return_value=(owner.pk, owner.version)):
            self.assertTrue(_node_ready_to_lease(nodes_by_id[owner.pk], parents_by_child, nodes_by_id))

    def test_no_cache_dependency_allows_node(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        nodes_by_id, parents_by_child = self._graph()
        self.assertTrue(_node_ready_to_lease(nodes_by_id[node.pk], parents_by_child, nodes_by_id))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestCascadeScheduling(AnalysisSetupMixin, TestCase):
    """A→B→C, all DIRTY: nodes are leased one layer at a time as parents settle.
    No wait_for_node signature is ever produced."""

    def _force_status(self, node, status):
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def setUp(self):
        self.a = AllVariantsNode.objects.create(analysis=self.analysis)
        self.b = GeneListNode.objects.create(analysis=self.analysis)
        self.c = GeneListNode.objects.create(analysis=self.analysis)
        self.b.add_parent(self.a)
        self.c.add_parent(self.b)

    def _dispatch(self):
        """Lease ready nodes and build their launch signatures (as create_and_launch would),
        asserting no wait_for_node task is ever produced."""
        leased = lease_ready_nodes(self.analysis.pk, "test-worker")
        for node_id, version in leased:
            for name in _signature_task_names(_node_launch_signature(node_id, version)):
                self.assertNotEqual(name, WAIT_FOR_NODE_TASK, "wait_for_node must never be scheduled")
        return leased

    def test_only_root_leased_first(self):
        self._dispatch()
        self.assertTrue(NodeTask.objects.filter(node_version__node=self.a).exists())
        self.assertFalse(NodeTask.objects.filter(node_version__node=self.b).exists())
        self.assertFalse(NodeTask.objects.filter(node_version__node=self.c).exists())

    def test_cascade_unblocks_one_layer_at_a_time(self):
        self._dispatch()  # leases A
        self._force_status(self.a, NodeStatus.READY)

        self._dispatch()  # A ready -> leases B (not C)
        self.assertTrue(NodeTask.objects.filter(node_version__node=self.b).exists())
        self.assertFalse(NodeTask.objects.filter(node_version__node=self.c).exists())
        self._force_status(self.b, NodeStatus.READY)

        self._dispatch()  # B ready -> leases C
        self.assertTrue(NodeTask.objects.filter(node_version__node=self.c).exists())


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDeadWorkerReclamation(AnalysisSetupMixin, TestCase):
    """A node whose owning worker died (expired lease, still QUEUED) is reclaimed by the
    single-worker dispatcher; after MAX_NODE_ATTEMPTS it is failed rather than hanging forever."""

    def _force_status(self, node, status):
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def _expire_lease(self, node):
        past = timezone.now() - timedelta(seconds=1)
        NodeTask.objects.filter(node_version__node=node).update(lease_expires=past)
        self._force_status(node, NodeStatus.QUEUED)

    def test_expired_lease_is_reclaimed_incrementing_attempts(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker-1")
        self.assertEqual(NodeTask.objects.get(node_version__node=node).attempt_count, 1)

        self._expire_lease(node)
        leased = lease_ready_nodes(self.analysis.pk, "worker-2")
        self.assertIn((node.pk, node.version), leased)
        self.assertEqual(NodeTask.objects.get(node_version__node=node).attempt_count, 2)

    def test_live_lease_is_not_reclaimed(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker-1")  # lease_expires is in the future
        leased = lease_ready_nodes(self.analysis.pk, "worker-2")
        self.assertEqual(leased, [])
        self.assertEqual(NodeTask.objects.get(node_version__node=node).attempt_count, 1)

    def test_terminal_fail_after_max_attempts(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        child = GeneListNode.objects.create(analysis=self.analysis)
        child.add_parent(node)

        for _ in range(MAX_NODE_ATTEMPTS):
            lease_ready_nodes(self.analysis.pk, "worker")
            self._expire_lease(node)

        # Next pass sees attempt_count == MAX_NODE_ATTEMPTS -> gives up (fails the parent)
        lease_ready_nodes(self.analysis.pk, "worker")
        node.refresh_from_db()
        self.assertEqual(node.status, NodeStatus.ERROR)
        # A subsequent dispatch now leases the child (parent settled to ERROR); it will itself
        # settle to ERROR_WITH_PARENT when it runs — no infinite hang.
        lease_ready_nodes(self.analysis.pk, "worker")
        self.assertTrue(NodeTask.objects.filter(node_version__node=child).exists())

    def test_reschedule_stalled_analyses_discovers_expired_lease(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker-1")
        self._expire_lease(node)
        with mock.patch.object(Signature, "apply_async") as m:
            reschedule_stalled_analyses()
        self.assertTrue(m.called, "stalled analysis should kick the dispatcher")


VENN_CACHE_COUNT_TASK = "analysis.models.nodes.filters.venn_node.venn_cache_count"


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestVennCacheScheduling(AnalysisSetupMixin, TestCase):
    """ VennNode builds its result in a VennNodeCache (separate from the regular NodeCache the
        dispatcher's gate understands). Under the issue #346 per-node launch model the cache count
        is chained ahead of the Venn's own update; these tests lock in that the chaining and the
        cache bookkeeping are idempotent across re-leases and shared parents, so a Venn update can
        never run against a half-built / CREATED cache. """

    def _force_status(self, node, status):
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def _make_venn(self, set_operation=SetOperations.INTERSECTION):
        a = AllVariantsNode.objects.create(analysis=self.analysis)
        b = AllVariantsNode.objects.create(analysis=self.analysis)
        venn = VennNode.objects.create(analysis=self.analysis, set_operation=set_operation)
        venn.add_parent(a, side=VennNode.LEFT_PARENT)
        venn.add_parent(b, side=VennNode.RIGHT_PARENT)
        venn.save()
        for n in (a, b):
            self._force_status(n, NodeStatus.READY)
        # Reload fresh so cached parents / ordered_parents reflect the committed edges + FKs
        venn = VennNode.objects.get(pk=venn.pk)
        self.assertTrue(venn.is_valid, f"Venn test node not valid: {venn.get_errors()}")
        return venn

    def _caches(self, venn):
        a, b = venn.ordered_parents
        return VennNodeCache.objects.filter(parent_a_node_version=a.node_version,
                                            parent_b_node_version=b.node_version)

    # --- launch signature: cache count is chained BEFORE the node's own update ---

    def test_launch_signature_chains_cache_before_update(self):
        venn = self._make_venn()
        names = _signature_task_names(_node_launch_signature(venn.pk, venn.version))
        self.assertIn(VENN_CACHE_COUNT_TASK, names)
        self.assertIn(venn.UPDATE_TASK, names)
        self.assertLess(names.index(VENN_CACHE_COUNT_TASK), names.index(venn.UPDATE_TASK),
                        "venn_cache_count must run before the Venn node's update")

    # --- idempotency of get_cache_task_args_set ---

    def test_relaunch_keeps_in_flight_cache(self):
        """ A re-lease (eg lease expiry) must NOT delete the CREATED collection an already-queued
            venn_cache_count is about to build - that destroy/recreate churn was what left a Venn
            update reading a 'Created' cache. """
        venn = self._make_venn()
        venn.get_cache_task_args_set()
        collection_pks = set(self._caches(venn).values_list("variant_collection_id", flat=True))
        self.assertTrue(all(collection_pks))

        task_args = venn.get_cache_task_args_set()  # simulate re-lease
        self.assertEqual(set(self._caches(venn).values_list("variant_collection_id", flat=True)),
                         collection_pks, "in-flight (CREATED) collection was replaced on re-lease")
        self.assertTrue(any(t for (t, _) in task_args),
                        "venn_cache_count must still be re-emitted while the cache is not SUCCESS")

    def test_success_cache_reused_emits_no_task(self):
        venn = self._make_venn()
        venn.get_cache_task_args_set()
        VariantCollection.objects.filter(
            pk__in=self._caches(venn).values_list("variant_collection_id", flat=True)
        ).update(status=ProcessingStatus.SUCCESS)

        task_args = venn.get_cache_task_args_set()
        self.assertTrue(all(t is None for (t, _) in task_args),
                        "A SUCCESS cache must be reused without re-emitting venn_cache_count")

    def test_error_cache_regenerated(self):
        venn = self._make_venn()
        venn.get_cache_task_args_set()
        old_pks = set(self._caches(venn).values_list("variant_collection_id", flat=True))
        VariantCollection.objects.filter(pk__in=old_pks).update(status=ProcessingStatus.ERROR)

        task_args = venn.get_cache_task_args_set()
        new_pks = set(self._caches(venn).values_list("variant_collection_id", flat=True))
        self.assertFalse(new_pks & old_pks, "ERROR collection should be discarded and replaced")
        self.assertTrue(any(t for (t, _) in task_args), "regenerated cache must re-emit venn_cache_count")

    # --- idempotency of the venn_cache_count task itself ---

    def test_venn_cache_count_skips_success(self):
        """ A sibling Venn sharing the same parents can schedule a second count for an
            already-built cache; the task must no-op rather than rebuild (which would duplicate
            records). Uses an obsolete parent node_version so any real build would blow up. """
        vc = VariantCollection.objects.create(name="venn-test", status=ProcessingStatus.SUCCESS)
        venn = self._make_venn()
        a, b = venn.ordered_parents
        cache = VennNodeCache.objects.create(parent_a_node_version=a.node_version,
                                             parent_b_node_version=b.node_version,
                                             intersection_type=VennNodeCache.INTERSECTION,
                                             variant_collection=vc)
        venn_cache_count(cache.pk)  # must return immediately, no rebuild
        vc.refresh_from_db()
        self.assertEqual(vc.status, ProcessingStatus.SUCCESS)

    def test_venn_cache_count_handles_missing_collection(self):
        venn = self._make_venn()
        a, b = venn.ordered_parents
        cache = VennNodeCache.objects.create(parent_a_node_version=a.node_version,
                                             parent_b_node_version=b.node_version,
                                             intersection_type=VennNodeCache.INTERSECTION,
                                             variant_collection=None)
        venn_cache_count(cache.pk)  # collection is None - must return without error

    def test_venn_cache_count_obsolete_cache_id(self):
        venn_cache_count(-1)  # missing VennNodeCache - must return without error

    def test_venn_cache_count_routed_to_analysis_workers(self):
        from django.conf import settings
        route = settings.CELERY_TASK_ROUTES.get(VENN_CACHE_COUNT_TASK)
        self.assertEqual(route, {"queue": "analysis_workers", "routing_key": "analysis_workers"},
                         "venn_cache_count must run on the analysis pool, not default db_workers")


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTransientErrorBackoff(AnalysisSetupMixin, TestCase):
    """Transient (OperationalError) failures back off via run_after and retry up to
    MAX_NODE_ATTEMPTS rather than dying immediately."""

    def _force_status(self, node, status):
        type(node).objects.filter(pk=node.pk).update(status=status)
        node.refresh_from_db()

    def test_next_backoff_sequence(self):
        self.assertEqual([next_backoff(0), next_backoff(1), next_backoff(2)], [10, 20, 40])

    def test_backoff_sets_dirty_and_future_run_after(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker")  # attempt_count -> 1, QUEUED
        with mock.patch.object(Signature, "apply_async"):
            self.assertTrue(_backoff_node(node.pk, node.version, self.analysis.pk))
        node.refresh_from_db()
        node_task = NodeTask.objects.get(node_version__node=node)
        self.assertEqual(node.status, NodeStatus.DIRTY)
        self.assertIsNotNone(node_task.run_after)
        self.assertGreater(node_task.run_after, timezone.now())
        self.assertIsNone(node_task.lease_expires)

    def test_node_in_backoff_window_is_not_released(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker")
        with mock.patch.object(Signature, "apply_async"):
            _backoff_node(node.pk, node.version, self.analysis.pk)  # DIRTY + future run_after
        leased = lease_ready_nodes(self.analysis.pk, "worker")
        self.assertEqual(leased, [], "node should not be leased while backing off")

    def test_node_released_after_backoff_window(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        lease_ready_nodes(self.analysis.pk, "worker")
        with mock.patch.object(Signature, "apply_async"):
            _backoff_node(node.pk, node.version, self.analysis.pk)
        NodeTask.objects.filter(node_version__node=node).update(run_after=timezone.now() - timedelta(seconds=1))
        leased = lease_ready_nodes(self.analysis.pk, "worker")
        self.assertIn((node.pk, node.version), leased)

    def test_backoff_gives_up_after_max_attempts(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        NodeTask.objects.create(node_version=node.node_version, analysis_update_uuid=uuid.uuid4(),
                                attempt_count=MAX_NODE_ATTEMPTS)
        with mock.patch.object(Signature, "apply_async"):
            self.assertFalse(_backoff_node(node.pk, node.version, self.analysis.pk))


class TestSingleWorkerInvariant(TestCase):
    """All node assignment/leasing happens only in lease_ready_nodes, reached only from
    create_and_launch_analysis_tasks. Worker tasks never call it."""

    @staticmethod
    def _task_source(task):
        # celery shared_task wraps the function; .run is the original callable
        return inspect.getsource(getattr(task, "run", task))

    def test_only_dispatcher_calls_lease_ready_nodes(self):
        # Check for an actual call site ("lease_ready_nodes("), not mentions in docstrings.
        self.assertIn("lease_ready_nodes(", self._task_source(create_and_launch_analysis_tasks))
        for task in (update_node_task, node_cache_task, reschedule_stalled_analyses):
            src = self._task_source(task)
            self.assertNotIn("lease_ready_nodes(", src,
                             f"{task.name} must not lease nodes — assignment is single-worker only")
