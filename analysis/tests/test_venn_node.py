from unittest.mock import patch

from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode, AnalysisNode
from analysis.models.enums import NodeStatus, SetOperations
from analysis.models.nodes.filters.venn_node import VennNode, VennNodeCache, venn_cache_count
from analysis.tests.utils import AnalysisSetupMixin
from snpdb.models import ProcessingStatus, VariantCollection

A = VennNodeCache.A_ONLY
I = VennNodeCache.INTERSECTION
B = VennNodeCache.B_ONLY


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestVennNode(AnalysisSetupMixin, TestCase):

    def _venn(self, set_operation=SetOperations.INTERSECTION):
        venn = VennNode.objects.create(analysis=self.analysis, set_operation=set_operation)
        parents = []
        for side in (VennNode.LEFT_PARENT, VennNode.RIGHT_PARENT):
            parent = AllVariantsNode.objects.create(analysis=self.analysis)
            AllVariantsNode.objects.filter(pk=parent.pk).update(status=NodeStatus.READY, count=1)
            parent.refresh_from_db()
            venn.add_parent(parent, side=side)
            parents.append(parent)
        venn.save()
        return venn, *parents

    @staticmethod
    def _reload(venn) -> VennNode:
        """ Fresh from DB, so left/right FKs and parent edges come off the committed rows """
        return VennNode.objects.get(pk=venn.pk)

    def test_set_operations_map_to_cache_sections(self):
        expected = {
            SetOperations.NONE: [],
            SetOperations.UNION: [A, I, B],
            SetOperations.A_NOT_B: [A],
            SetOperations.INTERSECTION: [I],
            SetOperations.A_ONLY: [A, I],
            SetOperations.B_NOT_A: [B],
            SetOperations.SYMMETRIC_DIFFERENCE: [A, B],
            SetOperations.B_ONLY: [I, B],
        }
        for set_operation, intersection_types in expected.items():
            node = VennNode(analysis=self.analysis, set_operation=set_operation)
            self.assertEqual(node.get_vennodecache_intersection_types(), intersection_types,
                             f"Cache sections for {set_operation.label}")

    def test_none_returns_nothing(self):
        venn, _left, _right = self._venn(set_operation=SetOperations.NONE)
        arg_q_dict = self._reload(venn)._get_node_cache_arg_q_dict()
        (q,) = arg_q_dict[None].values()
        self.assertEqual(q, venn.q_none())

    def test_none_has_no_cache_tasks_and_no_errors(self):
        """ Deselecting every region in the venn widget is a valid choice, not a broken node """
        venn, _left, _right = self._venn(set_operation=SetOperations.NONE)
        venn = self._reload(venn)
        self.assertEqual(venn.get_cache_task_args_set(), set())
        venn.refresh_from_db()
        self.assertNotEqual(venn.status, NodeStatus.ERROR)
        self.assertFalse(venn.errors)

    def test_parents_keep_their_sides(self):
        venn, left, right = self._venn()
        self.assertEqual(venn.get_side_for_parent(left), VennNode.LEFT_PARENT)
        self.assertEqual(venn.get_side_for_parent(right), VennNode.RIGHT_PARENT)

    def test_removing_a_parent_clears_its_side(self):
        venn, left, _right = self._venn()
        venn.remove_parent(left)
        self.assertIsNone(venn.left_parent)

    def test_side_lookup_survives_a_round_trip_through_the_db(self):
        """ The node editor loads both the venn and its parent fresh, as subclasses """
        venn, left, _right = self._venn()
        venn = self._reload(venn)
        parent = AnalysisNode.objects.get_subclass(pk=left.pk)
        self.assertEqual(venn.get_side_for_parent(parent), VennNode.LEFT_PARENT)

    def test_removing_a_parent_loaded_from_the_db_clears_its_side(self):
        venn, left, _right = self._venn()
        venn = self._reload(venn)
        venn.remove_parent(AnalysisNode.objects.get_subclass(pk=left.pk))
        venn.save()
        venn.refresh_from_db()
        self.assertIsNone(venn.left_parent_id)

    def test_ordered_parents_returns_left_then_right(self):
        venn, left, right = self._venn()
        a, b = self._reload(venn).ordered_parents
        self.assertEqual((a.pk, b.pk), (left.pk, right.pk))

    def _cache(self, venn, intersection_type) -> VennNodeCache:
        a, b = self._reload(venn).ordered_parents
        variant_collection = VariantCollection.objects.create(name="test venn cache")
        return VennNodeCache.objects.create(parent_a_node_version=a.node_version,
                                            parent_b_node_version=b.node_version,
                                            intersection_type=intersection_type,
                                            variant_collection=variant_collection)

    def test_intersection_with_an_empty_parent_skips_loading_variants(self):
        """ Pulling a side's variant ids is the expensive part - an empty parent decides it without them """
        venn, _left, right = self._venn(set_operation=SetOperations.INTERSECTION)
        AllVariantsNode.objects.filter(pk=right.pk).update(count=0)
        vennode_cache = self._cache(venn, VennNodeCache.INTERSECTION)

        with patch("analysis.models.nodes.filters.venn_node._node_variant_ids") as mock_variant_ids:
            venn_cache_count(vennode_cache.pk)
        mock_variant_ids.assert_not_called()

        vennode_cache.variant_collection.refresh_from_db()
        self.assertEqual(vennode_cache.variant_collection.status, ProcessingStatus.SUCCESS)

    def test_difference_only_skips_the_side_it_subtracts(self):
        """ A - B with an empty B is just A, so B's variants are never needed """
        venn, left, right = self._venn(set_operation=SetOperations.A_NOT_B)
        AllVariantsNode.objects.filter(pk=right.pk).update(count=0)
        vennode_cache = self._cache(venn, VennNodeCache.A_ONLY)

        with patch("analysis.models.nodes.filters.venn_node._node_variant_ids") as mock_variant_ids:
            mock_variant_ids.return_value = set()
            venn_cache_count(vennode_cache.pk)
        (call_node,) = [c.args[0] for c in mock_variant_ids.call_args_list]
        self.assertEqual(call_node.pk, left.pk)
