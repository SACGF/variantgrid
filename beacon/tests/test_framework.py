from django.test import TestCase, override_settings
from django.urls import reverse

from beacon.schema import OBSERVATIONS_DATASET_ID, CLASSIFICATIONS_DATASET_ID


@override_settings(BEACON_ENABLED=True)
class BeaconFrameworkTestCase(TestCase):
    """ Framework endpoints return spec-shaped JSON built from BEACON_CONFIG (§6, §11). """

    def test_info(self):
        response = self.client.get(reverse("beacon_info"))
        self.assertEqual(response.status_code, 200)
        payload = response.json()["response"]
        self.assertEqual(payload["apiVersion"], "v2.0.0")
        dataset_ids = {d["id"] for d in payload["datasets"]}
        self.assertEqual(dataset_ids, {OBSERVATIONS_DATASET_ID, CLASSIFICATIONS_DATASET_ID})

    def test_service_info(self):
        response = self.client.get(reverse("beacon_service_info"))
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(data["type"]["artifact"], "beacon")

    def test_configuration(self):
        response = self.client.get(reverse("beacon_configuration"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("entryTypes", response.json()["response"])

    def test_entry_types(self):
        response = self.client.get(reverse("beacon_entry_types"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("genomicVariant", response.json()["response"]["entryTypes"])

    def test_filtering_terms(self):
        response = self.client.get(reverse("beacon_filtering_terms"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("filteringTerms", response.json()["response"])

    def test_map(self):
        response = self.client.get(reverse("beacon_map"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("endpointSets", response.json()["response"])

    def test_disabled_returns_404(self):
        with override_settings(BEACON_ENABLED=False):
            response = self.client.get(reverse("beacon_info"))
        self.assertEqual(response.status_code, 404)
