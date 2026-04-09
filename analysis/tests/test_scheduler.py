"""
Tests for analysis node scheduling logic.

Covers NodeStatus categorisation and the core _get_analysis_update_tasks()
scheduler function. The expectedFailure test documents the broken behaviour
fixed by issue #346: remove @unittest.expectedFailure once the implementation
is done and verify it passes.
"""
import unittest

from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode, NodeTask
from analysis.models.enums import NodeStatus
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.tasks.analysis_update_tasks import _get_analysis_update_tasks
from analysis.tests.utils import AnalysisSetupMixin


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
    """Integration tests for _get_analysis_update_tasks().

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

    def test_no_nodes_returns_empty_task_list(self):
        tasks = _get_analysis_update_tasks(self.analysis.pk)
        self.assertEqual(tasks, [])

    def test_no_dirty_nodes_returns_empty_task_list(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        self._force_status(node, NodeStatus.READY)
        tasks = _get_analysis_update_tasks(self.analysis.pk)
        self.assertEqual(tasks, [])

    def test_single_dirty_node_returns_one_task(self):
        AllVariantsNode.objects.create(analysis=self.analysis)  # default status=DIRTY
        tasks = _get_analysis_update_tasks(self.analysis.pk)
        self.assertEqual(len(tasks), 1)

    def test_single_dirty_node_becomes_queued(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        _get_analysis_update_tasks(self.analysis.pk)
        node.refresh_from_db()
        self.assertEqual(node.status, NodeStatus.QUEUED)

    def test_single_dirty_node_gets_node_task_lock(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertTrue(NodeTask.objects.filter(node=node).exists())

    # ------------------------------------------------------------------
    # Double-dispatch: the NodeTask lock prevents scheduling the same node twice
    # ------------------------------------------------------------------

    def test_double_dispatch_creates_only_one_node_task(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        _get_analysis_update_tasks(self.analysis.pk)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertEqual(NodeTask.objects.filter(node=node).count(), 1)

    def test_double_dispatch_second_call_returns_no_tasks(self):
        AllVariantsNode.objects.create(analysis=self.analysis)
        _get_analysis_update_tasks(self.analysis.pk)
        tasks = _get_analysis_update_tasks(self.analysis.pk)
        self.assertEqual(tasks, [])

    # ------------------------------------------------------------------
    # DAG: parent-child scheduling
    # ------------------------------------------------------------------

    def test_ready_parent_allows_child_to_be_scheduled(self):
        """Child whose only parent is READY should be scheduled."""
        _, child = self._make_chain(parent_status=NodeStatus.READY)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertTrue(NodeTask.objects.filter(node=child).exists())

    def test_error_parent_allows_child_to_be_scheduled(self):
        """Error states are READY_STATUSES; child should still be scheduled
        (it will fail with ERROR_WITH_PARENT when it actually runs)."""
        _, child = self._make_chain(parent_status=NodeStatus.ERROR)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertTrue(NodeTask.objects.filter(node=child).exists())

    @unittest.expectedFailure  # Issue #346 — remove decorator once fix is implemented
    def test_loading_parent_child_not_locked(self):
        """Child should NOT be scheduled (locked via NodeTask) when its parent
        is still loading.  After the fix, the scheduler skips it and the
        post-completion trigger picks it up once the parent reaches READY.

        Currently FAILS: the scheduler locks the child and inserts wait_for_node
        into the chain, which is the worker-starvation pattern we are removing.
        """
        _, child = self._make_chain(parent_status=NodeStatus.LOADING)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertFalse(
            NodeTask.objects.filter(node=child).exists(),
            "Child was locked (NodeTask created) despite parent being LOADING — "
            "wait_for_node pattern is still active",
        )

    @unittest.expectedFailure  # Issue #346 — remove decorator once fix is implemented
    def test_queued_parent_child_not_locked(self):
        """QUEUED is a LOADING_STATUS: child must not be scheduled.
        This is the exact scenario from the VG.com crash log."""
        _, child = self._make_chain(parent_status=NodeStatus.QUEUED)
        _get_analysis_update_tasks(self.analysis.pk)
        self.assertFalse(
            NodeTask.objects.filter(node=child).exists(),
            "Child was locked despite parent being QUEUED",
        )
