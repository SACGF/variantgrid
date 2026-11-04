""" The analysis page inlines node data in a <script> block, so user-named nodes must not be able to end it """
import json

from django.template import Context, Template
from django.test import TestCase

from analysis.models import MergeNode
from analysis.tests.utils import AnalysisSetupMixin


class TestRenderNodesArray(AnalysisSetupMixin, TestCase):

    def test_node_name_cannot_close_script(self):
        name = "agent_test</script><b>x</b>"
        node = MergeNode.objects.create(analysis=self.analysis, name=name, auto_node_name=False)
        template = Template("{% load analysis_node_tags %}{% render_nodes_array nodes %}")
        output = template.render(Context({"nodes": [node]}))
        self.assertNotIn("<", output)
        self.assertEqual(name, json.loads(output)[0]["name"])
