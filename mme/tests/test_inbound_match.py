from types import SimpleNamespace
from unittest.mock import patch

from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.test import APIClient

from classification.enums.classification_enums import ClinicalSignificance
from classification.tests.models.test_utils import ClassificationTestUtils
from mme.matching import find_matches
from mme.models import MMEInboundMatch, MMEInboundQuery
from mme.tests.fakes import FakeClassification, make_classification

TOKEN = "inbound-secret"
QUERY_PATIENT = {"id": "q1", "features": [{"id": "HP:0001250"}],
                 "genomicFeatures": [{"gene": {"id": "BRCA1"}}]}


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN},
                   MME_DISCLAIMER="ours", MME_TERMS="our terms",
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class InboundMatchViewTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()
        self.client = APIClient()
        self.url = reverse("mme_api_match")

    def _post(self, body, token=TOKEN):
        headers = {}
        if token is not None:
            headers["HTTP_X_AUTH_TOKEN"] = token
        return self.client.post(self.url, data=body, format="json", **headers)

    def test_missing_token_rejected(self):
        response = self._post({"patient": QUERY_PATIENT}, token=None)
        self.assertEqual(response.status_code, 401)

    def test_wrong_token_rejected(self):
        response = self._post({"patient": QUERY_PATIENT}, token="nope")
        self.assertEqual(response.status_code, 401)

    def test_missing_features_and_genomic_rejected(self):
        response = self._post({"patient": {"id": "q1"}})
        self.assertEqual(response.status_code, 400)

    def test_valid_query_returns_results_and_writes_audit_row(self):
        response = self._post({"patient": QUERY_PATIENT})
        self.assertEqual(response.status_code, 200)
        self.assertIn("results", response.json())
        self.assertEqual(MMEInboundQuery.objects.count(), 1)
        self.assertEqual(MMEInboundQuery.objects.first().num_results, len(response.json()["results"]))

    def test_query_attributed_to_the_peer_that_authenticated(self):
        self._post({"patient": QUERY_PATIENT})
        self.assertEqual(MMEInboundQuery.objects.first().peer_node_id, "testpeer")

    def test_response_carries_our_disclaimer_and_terms(self):
        data = self._post({"patient": QUERY_PATIENT}).json()
        self.assertEqual(data["disclaimer"], "ours")
        self.assertEqual(data["terms"], "our terms")

    def test_returned_classifications_are_recorded_one_row_each(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        # An ineligible (Benign) record with the same gene must not be returned or recorded.
        make_classification(self.lab, self.user, gene_symbol="BRCA1",
                            clinical_significance=ClinicalSignificance.BENIGN)

        response = self._post({"patient": QUERY_PATIENT})
        self.assertEqual(len(response.json()["results"]), 2)

        inbound_query = MMEInboundQuery.objects.get()
        matches = MMEInboundMatch.objects.filter(inbound_query=inbound_query)
        self.assertEqual(matches.count(), 2)
        for match in matches:
            self.assertEqual(match.remote_patient_id, "q1")
            self.assertEqual(match.query_patient_json, QUERY_PATIENT)

    def test_test_mode_query_records_but_does_not_notify(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        with patch("mme.views_rest.notify_mme_depositors_task") as mock_task:
            self._post({"patient": dict(QUERY_PATIENT, test=True)})
        mock_task.si.assert_not_called()
        self.assertEqual(MMEInboundMatch.objects.count(), 1)

    def test_real_query_queues_depositor_notification(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        with patch("mme.views_rest.notify_mme_depositors_task") as mock_task:
            self._post({"patient": QUERY_PATIENT})
        mock_task.si.assert_called_once_with(MMEInboundQuery.objects.get().pk)

    @override_settings(MME_ENABLED=False)
    def test_disabled_rejected(self):
        response = self._post({"patient": QUERY_PATIENT})
        self.assertEqual(response.status_code, 401)


@override_settings(MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"},
                   MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
class FindMatchesTestCase(TestCase):
    """ Candidates come only from mme_eligible_classifications(), returned as MMEMatch so
        the view can both render and persist them. """

    def test_gene_overlap_scores_and_returns_public_classification(self):
        classification = FakeClassification(gene_symbol="BRCA1")
        fake_qs = SimpleNamespace(iterator=lambda: iter([SimpleNamespace(classification=classification)]))
        with patch("mme.matching.mme_eligible_classifications", return_value=fake_qs):
            matches = find_matches({"genomicFeatures": [{"gene": {"id": "BRCA1"}}]})

        self.assertEqual(len(matches), 1)
        self.assertGreater(matches[0].score, 0)
        self.assertIs(matches[0].classification, classification)
        result = matches[0].as_result()
        self.assertGreater(result["score"]["patient"], 0)
        self.assertEqual(result["patient"]["genomicFeatures"][0]["gene"]["id"], "BRCA1")

    def test_no_overlap_returns_nothing(self):
        fake_cm = SimpleNamespace(classification=FakeClassification(gene_symbol="TP53"))
        fake_qs = SimpleNamespace(iterator=lambda: iter([fake_cm]))
        with patch("mme.matching.mme_eligible_classifications", return_value=fake_qs):
            matches = find_matches({"genomicFeatures": [{"gene": {"id": "BRCA1"}}]})

        self.assertEqual(matches, [])
