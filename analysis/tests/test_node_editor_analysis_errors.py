from django.contrib.auth.models import User
from django.test import TestCase
from django.test.client import Client
from django.urls.base import reverse

from analysis.models import AllVariantsNode, Analysis
from analysis.tests.utils import AnalysisSetupMixin
from library.guardian_utils import assign_permission_to_user_and_groups


class NodeEditorAnalysisErrorsTest(AnalysisSetupMixin, TestCase):
    """ An analysis with broken settings (eg no annotation version) can't build node forms or grid columns,
        so the editor needs to load and report why rather than returning a 500 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get(username=f"test_{cls.__name__}")
        assign_permission_to_user_and_groups(cls.user, cls.analysis)
        cls.node = AllVariantsNode.objects.create(analysis=cls.analysis)

    def test_editor_reports_errors_when_annotation_version_missing(self):
        Analysis.objects.filter(pk=self.analysis.pk).update(annotation_version=None)

        client = Client()
        client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                           "analysis_version": self.analysis.version,
                                           "node_id": self.node.pk,
                                           "node_version": self.node.version,
                                           "extra_filters": "default"})
        response = client.get(url)
        self.assertEqual(200, response.status_code)
        errors = response.context["errors"]
        annotation_errors = [e for e in errors if "annotation_version" in e]
        self.assertTrue(annotation_errors, errors)
        for error in annotation_errors:
            self.assertNotIn("javascript:", error)  # Links belong in the template, not the message
