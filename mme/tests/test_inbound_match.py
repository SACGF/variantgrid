from types import SimpleNamespace
from unittest.mock import patch

from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.test import APIClient

from mme.matching import find_matches
from mme.models import MMEInboundQuery
from mme.tests.fakes import FakeClassification

TOKEN = "inbound-secret"
QUERY_PATIENT = {"id": "q1", "features": [{"id": "HP:0001250"}],
                 "genomicFeatures": [{"gene": {"id": "BRCA1"}}]}


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKEN=TOKEN,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class InboundMatchViewTestCase(TestCase):

    def setUp(self):
        self.client = APIClient()
        self.url = reverse("mme_api_match")

    def _post(self, body, token=TOKEN):
        headers = {}
        if token is not None:
            headers["HTTP_X_AUTH_TOKEN"] = token
        return self.client.post(self.url, data=body, format="json", **headers)

    def test_missing_token_rejected(self):
        response = self._post({"patient": QUERY_PATIENT}, token=None)
        self.assertIn(response.status_code, (401, 403))

    def test_wrong_token_rejected(self):
        response = self._post({"patient": QUERY_PATIENT}, token="nope")
        self.assertIn(response.status_code, (401, 403))

    def test_missing_features_and_genomic_rejected(self):
        response = self._post({"patient": {"id": "q1"}})
        self.assertEqual(response.status_code, 400)

    def test_valid_query_returns_results_and_writes_audit_row(self):
        response = self._post({"patient": QUERY_PATIENT})
        self.assertEqual(response.status_code, 200)
        self.assertIn("results", response.json())
        self.assertEqual(MMEInboundQuery.objects.count(), 1)
        self.assertEqual(MMEInboundQuery.objects.first().num_results, len(response.json()["results"]))

    @override_settings(MME_ENABLED=False)
    def test_disabled_rejected(self):
        response = self._post({"patient": QUERY_PATIENT})
        self.assertIn(response.status_code, (401, 403))


@override_settings(MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"},
                   MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
class FindMatchesTestCase(TestCase):
    """ find_matches sources candidates only from mme_eligible_classifications() (PUBLIC). """

    def test_gene_overlap_scores_and_returns_public_classification(self):
        fake_cm = SimpleNamespace(classification=FakeClassification(gene_symbol="BRCA1"))
        fake_qs = SimpleNamespace(iterator=lambda: iter([fake_cm]))
        with patch("mme.matching.mme_eligible_classifications", return_value=fake_qs):
            results = find_matches({"genomicFeatures": [{"gene": {"id": "BRCA1"}}]})

        self.assertEqual(len(results), 1)
        self.assertGreater(results[0]["score"]["patient"], 0)
        self.assertEqual(results[0]["patient"]["genomicFeatures"][0]["gene"]["id"], "BRCA1")

    def test_no_overlap_returns_nothing(self):
        fake_cm = SimpleNamespace(classification=FakeClassification(gene_symbol="TP53"))
        fake_qs = SimpleNamespace(iterator=lambda: iter([fake_cm]))
        with patch("mme.matching.mme_eligible_classifications", return_value=fake_qs):
            results = find_matches({"genomicFeatures": [{"gene": {"id": "BRCA1"}}]})

        self.assertEqual(results, [])
