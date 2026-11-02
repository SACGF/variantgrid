"""
AnalysisNode.node_type and the select_subclasses() it narrows: only the subclass tables of the node types a queryset
returns are joined (#2586).
"""
from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from analysis.models import AnalysisNode, MergeNode, VennNode
from analysis.models.nodes.node_types import get_node_types_hash
from analysis.tests.utils import AnalysisSetupMixin


class NodeTypeTest(AnalysisSetupMixin, TestCase):

    def test_every_node_class_loads_as_itself(self):
        nodes = [node_class.objects.create(analysis=self.analysis) for node_class in get_node_types_hash().values()]
        clone = nodes[0].save_clone()
        expected = {node.pk: type(node) for node in [*nodes, clone]}

        loaded = {node.pk: type(node) for node in AnalysisNode.objects.filter(pk__in=expected).select_subclasses()}
        self.assertEqual(expected, loaded)
        self.assertEqual({cls.__name__ for cls in expected.values()},
                         set(AnalysisNode.objects.filter(pk__in=expected).values_list("node_type", flat=True)))

    def test_joins_only_the_present_subclass_tables(self):
        merge = MergeNode.objects.create(analysis=self.analysis)
        venn = VennNode.objects.create(analysis=self.analysis)

        with CaptureQueriesContext(connection) as ctx:
            loaded = list(AnalysisNode.objects.filter(pk__in=[merge.pk, venn.pk]).order_by("pk").select_subclasses())
        self.assertEqual([MergeNode, VennNode], [type(node) for node in loaded])
        node_types_sql, nodes_sql = (q["sql"] for q in ctx.captured_queries)
        self.assertIn("DISTINCT", node_types_sql)
        self.assertIn('"analysis_mergenode"', nodes_sql)
        self.assertIn('"analysis_vennnode"', nodes_sql)
        self.assertNotIn('"analysis_samplenode"', nodes_sql)

    def test_get_slice_and_empty(self):
        venn = VennNode.objects.create(analysis=self.analysis)
        base = AnalysisNode.objects.get(pk=venn.pk)

        with CaptureQueriesContext(connection) as ctx:
            self.assertIsInstance(base.get_subclass(), VennNode)  # node_type already known
        self.assertEqual(1, len(ctx.captured_queries))

        self.assertIsInstance(AnalysisNode.objects.get_subclass(pk=venn.pk), VennNode)
        self.assertIsInstance(VennNode.objects.get_subclass(pk=venn.pk), VennNode)  # a subclass's own manager
        self.assertIsInstance(AnalysisNode.objects.filter(pk=venn.pk).select_subclasses().first(), VennNode)

        with CaptureQueriesContext(connection) as ctx:
            self.assertEqual([], list(AnalysisNode.objects.filter(pk=-1).select_subclasses()))
        self.assertEqual(1, len(ctx.captured_queries))
