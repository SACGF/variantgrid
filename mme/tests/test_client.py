from unittest.mock import patch, MagicMock

import requests
from django.test import TestCase, override_settings

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from mme.client import submit
from mme.models import MMESubmission, MMESubmissionStatus, MMEMatchResult

NODES = {"testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1"}}


@override_settings(MME_NODES=NODES, MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"},
                   MME_DISCLAIMER="please confirm", MME_ENABLED_PRODUCTION_SUBMIT=False)
class ClientTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, user = ClassificationTestUtils.lab_and_user()
        self.classification = Classification.create(
            user=user, lab=lab, data={SpecialEKeys.GENE_SYMBOL: {'value': 'BRCA1'}},
            save=True, source=SubmissionSource.API)
        self.submission = MMESubmission.objects.create(
            classification=self.classification, node_id="testnode", external_patient_id="vg:1")

    def _mock_response(self, payload):
        response = MagicMock()
        response.raise_for_status.return_value = None
        response.json.return_value = payload
        return response

    def test_success_sends_headers_and_persists_results(self):
        payload = {"results": [
            {"score": {"patient": 0.83},
             "patient": {"id": "remote-1", "contact": {"name": "Lab X", "href": "mailto:x@x.org"}}},
        ]}
        with patch("mme.client.requests.post", return_value=self._mock_response(payload)) as mock_post:
            submit(self.submission)

        _args, kwargs = mock_post.call_args
        headers = kwargs["headers"]
        self.assertEqual(headers["X-Auth-Token"], "secret-token")
        self.assertEqual(headers["Content-Type"], "application/vnd.ga4gh.matchmaker.v1.1+json")
        self.assertEqual(headers["Accept"], "application/vnd.ga4gh.matchmaker.v1.1+json")
        self.assertEqual(kwargs["url"], "https://node.test/match")
        self.assertEqual(kwargs["json"]["disclaimer"], "please confirm")

        self.submission.refresh_from_db()
        self.assertEqual(self.submission.status, MMESubmissionStatus.SUBMITTED)
        self.assertIsNotNone(self.submission.submitted)
        result = MMEMatchResult.objects.get(submission=self.submission)
        self.assertEqual(result.score, 0.83)
        self.assertEqual(result.matched_patient_id, "remote-1")
        self.assertEqual(result.contact_name, "Lab X")

    def test_failure_sets_error_and_notifies_admin(self):
        response = MagicMock()
        response.raise_for_status.side_effect = requests.HTTPError("500 Server Error")
        with patch("mme.client.requests.post", return_value=response), \
                patch("mme.client.AdminNotificationBuilder") as mock_nb:
            with self.assertRaises(requests.HTTPError):
                submit(self.submission)

        self.submission.refresh_from_db()
        self.assertEqual(self.submission.status, MMESubmissionStatus.ERROR)
        self.assertIn("500 Server Error", self.submission.error)
        mock_nb.assert_called_once()
        mock_nb.return_value.send.assert_called_once()
        self.assertFalse(MMEMatchResult.objects.filter(submission=self.submission).exists())
