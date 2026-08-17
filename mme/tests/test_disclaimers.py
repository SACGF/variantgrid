from unittest.mock import MagicMock, patch

from django.test import TestCase, override_settings
from django.urls import reverse

from classification.tests.models.test_utils import ClassificationTestUtils
from library.django_utils.unittest_utils import URLTestCase
from mme.client import submit
from mme.disclaimers import (
    connected_nodes,
    effective_disclaimer,
    mme_response_body,
    node_disclaimer,
)
from mme.models import MMESubmission
from mme.tests.fakes import make_classification

NODES = {
    "testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1",
                 "disclaimer": "baseline disclaimer", "terms": "baseline terms"},
    "quietnode": {"base_url": "https://quiet.test", "token": "t", "api_version": "1.1"},
}


@override_settings(MME_NODES=NODES)
class DisclaimerResolutionTestCase(TestCase):

    def test_baseline_comes_from_node_settings(self):
        self.assertEqual(node_disclaimer("testnode"),
                         {"disclaimer": "baseline disclaimer", "terms": "baseline terms"})

    def test_node_without_disclaimer_or_unknown_node_yields_blanks(self):
        self.assertEqual(node_disclaimer("quietnode"), {"disclaimer": "", "terms": ""})
        self.assertEqual(node_disclaimer("nosuchnode"), {"disclaimer": "", "terms": ""})

    def test_live_response_supersedes_baseline(self):
        effective = effective_disclaimer("testnode", "live disclaimer", "live terms")
        self.assertEqual(effective["disclaimer"], "live disclaimer")
        self.assertEqual(effective["terms"], "live terms")
        self.assertEqual(effective["source"], "response")

    def test_falls_back_to_baseline_when_response_carries_neither(self):
        effective = effective_disclaimer("testnode", "", "")
        self.assertEqual(effective["disclaimer"], "baseline disclaimer")
        self.assertEqual(effective["terms"], "baseline terms")
        self.assertEqual(effective["source"], "baseline")

    def test_connected_nodes_lists_every_configured_node(self):
        self.assertEqual([n["node_id"] for n in connected_nodes()], ["quietnode", "testnode"])

    @override_settings(MME_DISCLAIMER="ours", MME_TERMS="our terms")
    def test_response_body_carries_our_disclaimer_and_terms(self):
        self.assertEqual(mme_response_body({"results": []}),
                         {"results": [], "disclaimer": "ours", "terms": "our terms"})

    @override_settings(MME_DISCLAIMER="", MME_TERMS="")
    def test_response_body_omits_unset_disclaimer_and_terms(self):
        self.assertEqual(mme_response_body({"results": []}), {"results": []})


@override_settings(MME_NODES=NODES, MME_ENABLED=True,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class SubmissionDisclaimerTestCase(URLTestCase):
    """ What a node returned with THIS response is persisted and shown with these results. """

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, user = ClassificationTestUtils.lab_and_user()
        lab.contact_email = "curator@lab.org"
        lab.mme_enabled = True
        lab.save()
        self.user = user
        self.classification = make_classification(lab, user)
        self.submission = MMESubmission.objects.create(
            classification=self.classification, node_id="testnode", external_patient_id="vg:1")

    def _submit_with_response(self, payload):
        response = MagicMock()
        response.raise_for_status.return_value = None
        response.json.return_value = payload
        with patch("mme.client.requests.post", return_value=response):
            submit(self.submission)
        self.submission.refresh_from_db()

    def test_submit_persists_returned_disclaimer_and_terms(self):
        self._submit_with_response({"results": [], "disclaimer": "live d", "terms": "live t"})
        self.assertEqual(self.submission.response_disclaimer, "live d")
        self.assertEqual(self.submission.response_terms, "live t")

    def test_submission_page_renders_the_live_disclaimer(self):
        self._submit_with_response({"results": [], "disclaimer": "live d", "terms": "live t"})
        self.client.force_login(self.user)
        response = self.client.get(reverse("mme_view_submission", kwargs={"submission_id": self.submission.pk}))
        self.assertContains(response, "live d")
        self.assertNotContains(response, "baseline disclaimer")

    def test_submission_page_falls_back_to_the_published_baseline(self):
        self._submit_with_response({"results": []})
        self.client.force_login(self.user)
        response = self.client.get(reverse("mme_view_submission", kwargs={"submission_id": self.submission.pk}))
        self.assertContains(response, "baseline disclaimer")
        self.assertContains(response, "no disclaimer returned with this response")


@override_settings(MME_NODES=NODES)
class PublicDisclaimersPageTestCase(URLTestCase):
    """ Posted on our website, readable without an account. """

    def test_page_lists_every_connected_node_anonymously(self):
        response = self.client.get(reverse("mme_public_disclaimers"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "baseline disclaimer")
        self.assertContains(response, "baseline terms")
        self.assertContains(response, "quietnode")
