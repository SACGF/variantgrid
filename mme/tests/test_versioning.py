""" MME carries its API version in the media type, and a matching MAJOR version is
compatible. Most of the federation speaks v1.0, so a v1.0 peer must be served by our v1.1
implementation - see mme/versioning.py. """
import json

from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.request import Request
from rest_framework.test import APIClient, APIRequestFactory

from classification.tests.models.test_utils import ClassificationTestUtils
from mme.versioning import MME_MEDIA_TYPE, compatible_media_type, request_version

TOKEN = "inbound-secret"
QUERY_PATIENT = {"id": "q1", "features": [{"id": "HP:0001250"}],
                 "genomicFeatures": [{"gene": {"id": "BRCA1"}}]}

V1_0 = "application/vnd.ga4gh.matchmaker.v1.0+json"
V1_1 = "application/vnd.ga4gh.matchmaker.v1.1+json"
V2_0 = "application/vnd.ga4gh.matchmaker.v2.0+json"


class CompatibleMediaTypeTestCase(TestCase):

    def test_same_major_version_is_compatible(self):
        for media_type in (V1_0, V1_1, "application/vnd.ga4gh.matchmaker.v1.7+json"):
            self.assertEqual(compatible_media_type(media_type), media_type)

    def test_charset_parameter_ignored(self):
        self.assertEqual(compatible_media_type(f"{V1_0}; charset=utf-8"), V1_0)

    def test_other_major_version_incompatible(self):
        self.assertIsNone(compatible_media_type(V2_0))

    def test_non_mme_media_types_incompatible(self):
        for media_type in ("application/json", "text/html", "", None):
            self.assertIsNone(compatible_media_type(media_type))


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN},
                   MME_DISCLAIMER="ours", MME_CONTACT={"name": "T", "href": "mailto:t@t.org"})
class InboundVersionNegotiationTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.client = APIClient()
        self.url = reverse("mme_api_match")

    def _post(self, content_type, accept=None):
        return self.client.post(
            self.url, data=json.dumps({"patient": QUERY_PATIENT}),
            content_type=content_type, HTTP_ACCEPT=accept or content_type,
            HTTP_X_AUTH_TOKEN=TOKEN)

    def test_v1_0_peer_is_served(self):
        """ Six of the eight deployed nodes are v1.0 only - rejecting them is not an option. """
        response = self._post(V1_0)
        self.assertEqual(response.status_code, 200)
        self.assertIn("results", response.json())

    def test_v1_0_response_echoes_requested_version(self):
        response = self._post(V1_0)
        self.assertEqual(response["Content-Type"].split(";")[0], V1_0)

    def test_v1_1_peer_is_served(self):
        response = self._post(V1_1)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"].split(";")[0], V1_1)

    def test_plain_json_still_accepted(self):
        response = self._post("application/json")
        self.assertEqual(response.status_code, 200)

    def test_unsupported_major_version_is_406_not_415(self):
        response = self._post(V2_0, accept="*/*")
        self.assertEqual(response.status_code, 406)

    def test_406_advertises_our_latest_supported_version(self):
        response = self._post(V2_0, accept="*/*")
        self.assertEqual(response["Content-Type"].split(";")[0], MME_MEDIA_TYPE)

    def test_error_body_has_message(self):
        """ Spec: error responses carry a human-readable `message` (DRF defaults to `detail`). """
        response = self._post(V2_0, accept="*/*")
        self.assertIn("message", response.json())


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN}, MME_DISCLAIMER="ours")
class MetricsVersionNegotiationTestCase(TestCase):
    """ /metrics lists its renderers plain-JSON first (metrics-api.md specifies
        application/json), so version negotiation is worth asserting separately. """

    def setUp(self):
        self.client = APIClient()
        self.url = reverse("mme_api_metrics")

    def _get(self, accept):
        return self.client.get(self.url, HTTP_ACCEPT=accept, HTTP_X_AUTH_TOKEN=TOKEN)

    def test_v1_0_accept_is_served_and_echoed(self):
        response = self._get(V1_0)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"].split(";")[0], V1_0)

    def test_default_accept_is_plain_json(self):
        response = self._get("*/*")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"].split(";")[0], "application/json")

    def test_unsupported_major_version_is_406(self):
        response = self._get(V2_0)
        self.assertEqual(response.status_code, 406)
        self.assertEqual(response["Content-Type"].split(";")[0], MME_MEDIA_TYPE)


class RequestVersionTestCase(TestCase):
    """ The seam a future major would branch on. """

    def _request(self, content_type=None, accept=None):
        extra = {}
        if accept:
            extra["HTTP_ACCEPT"] = accept
        return Request(APIRequestFactory().post(
            "/mme/api/match", data="{}", content_type=content_type or "application/json", **extra))

    def test_reads_major_minor_from_content_type(self):
        self.assertEqual(request_version(self._request(content_type=V1_0)), ("1", "0"))

    def test_falls_back_to_accept_header(self):
        request = self._request(content_type="application/json", accept=V1_1)
        self.assertEqual(request_version(request), ("1", "1"))

    def test_plain_json_defaults_to_the_version_we_implement(self):
        self.assertEqual(request_version(self._request()), ("1", "1"))
