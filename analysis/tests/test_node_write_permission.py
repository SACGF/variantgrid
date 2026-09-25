from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse
from guardian.shortcuts import assign_perm

from analysis.models import AllVariantsNode, AnalysisNode
from analysis.tests.utils import AnalysisSetupMixin
from library.guardian_utils import DjangoPermission


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestNodeWritePermission(AnalysisSetupMixin, TestCase):
    """ Views that add nodes or cancel a load change the analysis - viewing it is not enough """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.viewer = User.objects.get_or_create(username="test_node_write_permission_viewer")[0]
        assign_perm(DjangoPermission.perm(cls.analysis, DjangoPermission.READ), cls.viewer, cls.analysis)
        cls.node = AllVariantsNode.objects.create(analysis=cls.analysis)

    def test_read_only_user_cannot_change_analysis(self):
        self.client.force_login(self.viewer)
        kwargs = {"analysis_id": self.analysis.pk, "node_id": self.node.pk}
        num_nodes = AnalysisNode.objects.filter(analysis=self.analysis).count()
        for url_name, data in [("create_filter_child", {"column_name": "variant", "column_filter": "x"}),
                               ("create_selected_child", {}),
                               ("node_cancel_load", {})]:
            with self.subTest(url_name=url_name):
                response = self.client.post(reverse(url_name, kwargs=kwargs), data)
                self.assertEqual(response.status_code, 403)
        self.assertEqual(AnalysisNode.objects.filter(analysis=self.analysis).count(), num_nodes)
