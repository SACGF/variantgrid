"""
A node's warnings and errors reach the analysis canvas through nodes_status, so the card can show them
without the editor being opened (#347).
"""
import json
from unittest.mock import patch

from django.test import Client
from django.urls import reverse

from analysis.models import NodeStatus
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.tests.test_grid_export import GridExportTestCase


class TestNodeStatusWarnings(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.client = Client()
        self.client.force_login(self.user)

    def _node_status(self, node) -> dict:
        url = reverse("nodes_status", kwargs={"analysis_id": self.analysis.pk})
        response = self.client.get(url, {"nodes": json.dumps([node.pk])})
        return response.json()["node_status"][0]

    def test_load_snapshots_warnings_for_the_card(self):
        node = self._sample_node()
        with patch.object(SampleNode, "get_warnings", return_value=["Genes of interest have incomplete coverage"]):
            node.load()

        self.assertEqual(["Genes of interest have incomplete coverage"], self._node_status(node)["warnings"])

    def test_errored_node_tooltip_hides_traceback(self):
        node = self._sample_node()
        traceback = 'Traceback (most recent call last):\n  File "x.py", line 1\nValueError: bad zygosity\n'
        node.update(status=NodeStatus.ERROR, errors=traceback)

        node_status = self._node_status(node)
        self.assertFalse(node_status["valid"])
        self.assertEqual(["Internal Error"], node_status["errors"])
