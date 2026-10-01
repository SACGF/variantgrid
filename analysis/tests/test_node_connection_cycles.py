import json
from unittest.mock import patch

from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.models.nodes.analysis_node import AnalysisDag, AnalysisEdge
from analysis.models.nodes.filters.filter_node import FilterNode
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.sources.all_variants_node import AllVariantsNode
from analysis.tests.utils import AnalysisSetupMixin


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
@patch("analysis.views.views_json.update_analysis")
class TestNodeConnectionCycles(AnalysisSetupMixin, TestCase):
    """ source -> venn (left) -> child -> grandchild, venn's right input free (#2060) """

    def setUp(self):
        self.client.force_login(self.analysis.user)
        self.source = AllVariantsNode.objects.create(analysis=self.analysis)
        self.venn = VennNode.objects.create(analysis=self.analysis)
        self.venn.add_parent(self.source, side=VennNode.LEFT_PARENT)
        self.venn.save()
        self.child = FilterNode.objects.create(analysis=self.analysis)
        self.child.add_parent(self.venn)
        self.grandchild = FilterNode.objects.create(analysis=self.analysis)
        self.grandchild.add_parent(self.child)

    def _connect(self, parent, node, **params):
        url = reverse("node_update", kwargs={"analysis_id": self.analysis.pk, "node_id": node.pk})
        params["parent_id"] = parent.pk
        response = self.client.post(url, {"op": "update_connection", "params": json.dumps(params)})
        self.assertEqual(response.status_code, 200)
        return response.json()

    def _edge_exists(self, parent, child) -> bool:
        return AnalysisEdge.objects.filter(parent=parent, child=child).exists()

    def test_cycle_rejected(self, _update_analysis):
        for node, params in [(self.venn, {"side": VennNode.RIGHT_PARENT}),
                             (self.child, {})]:
            with self.subTest(node=node.__class__.__name__):
                data = self._connect(self.grandchild, node, **params)
                self.assertTrue(data.get("non_fatal"))
                self.assertFalse(self._edge_exists(self.grandchild, node))

        self.assertIsNone(VennNode.objects.get(pk=self.venn.pk).right_parent)

    def test_diamond_allowed(self, _update_analysis):
        """ Two paths up to the same ancestor is not a cycle """
        merge = MergeNode.objects.create(analysis=self.analysis)
        merge.add_parent(self.grandchild)

        data = self._connect(self.child, merge)
        self.assertNotIn("non_fatal", data)
        self.assertTrue(self._edge_exists(self.child, merge))

    def test_dag_walks(self, _update_analysis):
        second_source = AllVariantsNode.objects.create(analysis=self.analysis)
        self.child.add_parent(second_source)
        dag = AnalysisDag(self.analysis.pk)

        self.assertEqual(dag.root_ids(self.grandchild.pk), {self.source.pk, second_source.pk})
        self.assertEqual(dag.root_ids(self.source.pk), set())
        self.assertEqual(dag.descendant_ids(self.venn.pk), {self.child.pk, self.grandchild.pk})
