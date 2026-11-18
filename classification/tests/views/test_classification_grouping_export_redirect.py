"""The grouping export summary page: what it says about the chosen settings, and how it handles bad ones."""
from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from classification.tests.models.test_utils import ClassificationTestUtils


class TestClassificationGroupingExportRedirect(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.client.force_login(self.user)

    def test_summary_and_refresh(self):
        query_string = "format=csv&labs=include-only-my-labs&rows_per_file=50&since=2026-01-01"
        with CaptureQueriesContext(connection) as ctx:
            response = self.client.get(reverse("classification_grouping_export_redirect") + "?" + query_string)
        self.assertEqual(response.status_code, 200)
        # describing the export must not run it
        self.assertFalse([q for q in ctx.captured_queries if "classification_classificationgrouping" in q["sql"]])

        summary = response.context["summary"]
        self.assertEqual(summary.included_labs, [self.lab])
        self.assertEqual(summary.rows_per_file, 100)
        self.assertTrue(summary.filename.endswith(".zip"))

        export_url = reverse("classification_grouping_export") + "?" + query_string
        self.assertContains(response, f'content="1;url={export_url.replace("&", "&amp;")}"')

    def test_invalid_format_shows_error(self):
        response = self.client.get(reverse("classification_grouping_export_redirect") + "?format=xls&labs=include-all-labs")
        self.assertEqual(response.status_code, 200)
        self.assertIsNone(response.context["summary"])
        self.assertIn("xls", response.context["error"])
        self.assertNotContains(response, 'http-equiv="refresh"')
