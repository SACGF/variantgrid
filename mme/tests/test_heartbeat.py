""" MME Heartbeat API - https://github.com/ga4gh/mme-apis/blob/1.1.1/heartbeat-api.md """
from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.test import APIClient

from manual.models import Deployment
from mme.versioning import supported_media_types

TOKEN = "inbound-secret"


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN},
                   MME_TEST=False, MME_DISCLAIMER="ours", MME_TERMS="our terms")
class HeartbeatEndpointTestCase(TestCase):

    def setUp(self):
        self.client = APIClient()
        self.url = reverse("mme_api_heartbeat")

    def _get(self, token=TOKEN):
        headers = {"HTTP_X_AUTH_TOKEN": token} if token is not None else {}
        return self.client.get(self.url, **headers)

    def test_requires_a_peer_token(self):
        self.assertEqual(self._get(token=None).status_code, 401)

    def test_wrong_token_rejected(self):
        self.assertEqual(self._get(token="nope").status_code, 401)

    @override_settings(MME_ENABLED=False, MME_INBOUND_TOKENS={})
    def test_answers_while_mme_is_disabled(self):
        """ Liveness is the point - a prospective peer checks reachability and supported
            versions before any token exists to authenticate with. """
        response = self._get(token=None)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["heartbeat"]["accept"], supported_media_types())

    def test_returns_mandatory_heartbeat_fields(self):
        response = self._get()
        self.assertEqual(response.status_code, 200)
        heartbeat = response.json()["heartbeat"]
        self.assertEqual(heartbeat["production"], True)
        self.assertIn("version", heartbeat)
        self.assertEqual(heartbeat["accept"], supported_media_types())

    def test_accept_advertises_both_versions_we_serve(self):
        """ The point of the endpoint: a peer negotiates from this rather than guessing. """
        accept = self._get().json()["heartbeat"]["accept"]
        self.assertIn("application/vnd.ga4gh.matchmaker.v1.0+json", accept)
        self.assertIn("application/vnd.ga4gh.matchmaker.v1.1+json", accept)

    @override_settings(MME_TEST=True)
    def test_test_instance_says_so(self):
        self.assertEqual(self._get().json()["heartbeat"]["production"], False)

    def test_carries_our_disclaimer_and_terms(self):
        data = self._get().json()
        self.assertEqual(data["disclaimer"], "ours")
        self.assertEqual(data["terms"], "our terms")

    def test_version_reports_last_deployment(self):
        Deployment.objects.create(git_hash="abcdef1234567890")
        self.assertEqual(self._get().json()["heartbeat"]["version"], "abcdef12")

    def test_version_without_a_recorded_deployment(self):
        self.assertEqual(self._get().json()["heartbeat"]["version"], "unknown")
