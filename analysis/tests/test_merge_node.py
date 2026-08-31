from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode
from analysis.models.enums import NodeStatus
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.tests.utils import AnalysisSetupMixin


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestMergeNode(AnalysisSetupMixin, TestCase):
    """ Merging parents that ask the same question is a no-op """

    def _ready_all_variants_node(self, **kwargs) -> AllVariantsNode:
        node = AllVariantsNode.objects.create(analysis=self.analysis, **kwargs)
        AllVariantsNode.objects.filter(pk=node.pk).update(status=NodeStatus.READY, count=1)
        node.refresh_from_db()
        return node

    def _merge_of(self, *parents) -> MergeNode:
        merge = MergeNode.objects.create(analysis=self.analysis)
        for parent in parents:
            merge.add_parent(parent)
        merge.save()
        return MergeNode.objects.get(pk=merge.pk)

    def test_identical_parents_do_not_modify(self):
        a = self._ready_all_variants_node()
        b = self._ready_all_variants_node()
        self.assertFalse(self._merge_of(a, b).modifies_parents())

    def test_different_parents_modify(self):
        a = self._ready_all_variants_node()
        b = self._ready_all_variants_node(max_het_or_hom_count=0)
        self.assertTrue(self._merge_of(a, b).modifies_parents())

    def test_identical_parents_collapse_to_a_single_parent(self):
        """ Lets the merge re-use a parent's grid/cache instead of running its own query """
        a = self._ready_all_variants_node()
        b = self._ready_all_variants_node()
        merge = self._merge_of(a, b)
        self.assertIn(merge.get_single_parent().pk, {a.pk, b.pk})
