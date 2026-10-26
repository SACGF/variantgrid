""" nodes_delete removes a selection of nodes bottom-up in one transaction (#2475) """
import json
from unittest.mock import patch

from django.test import TestCase
from django.urls import reverse

from analysis.models import AnalysisNode, FilterNode, TrioNode
from analysis.tests.utils import AnalysisSetupMixin


class TestNodesDelete(AnalysisSetupMixin, TestCase):

    def setUp(self):
        self.parent = TrioNode.objects.create(analysis=self.analysis)
        self.child = FilterNode.objects.create(analysis=self.analysis)
        self.grandchild = FilterNode.objects.create(analysis=self.analysis)
        self.parent.add_child(self.child)
        self.child.add_child(self.grandchild)
        self.client.force_login(self.analysis.user)

    def _delete(self, *nodes):
        url = reverse('nodes_delete', kwargs={"analysis_id": self.analysis.pk})
        return self.client.post(url, {"nodes": json.dumps([n.pk for n in nodes])})

    def test_failure_part_way_rolls_back(self):
        # Bottom-up, so the child is deleted and the grandchild detached before the parent's delete fails
        with patch.object(TrioNode, "delete", side_effect=RuntimeError("delete failed")):
            with self.assertRaises(RuntimeError):
                self._delete(self.parent, self.child)

        self.assertEqual(3, AnalysisNode.objects.filter(analysis=self.analysis).count())
        grandchild = FilterNode.objects.get(pk=self.grandchild.pk)
        self.assertEqual([self.child.pk], [p.pk for p in grandchild.parents()])

    def test_delete_detaches_kept_child(self):
        with patch("analysis.views.views_json.update_analysis") as update:
            with self.captureOnCommitCallbacks(execute=True):
                self.assertEqual(200, self._delete(self.parent, self.child).status_code)
                update.assert_not_called()
            update.assert_called_once_with(self.analysis.pk)

        self.assertEqual([self.grandchild.pk], list(AnalysisNode.objects.filter(analysis=self.analysis)
                                                    .values_list("pk", flat=True)))
        self.assertEqual([], list(FilterNode.objects.get(pk=self.grandchild.pk).parents()))
