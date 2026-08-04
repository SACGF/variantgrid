from django.core.cache import cache
from django.test import TestCase, override_settings
from django.urls import reverse
from rest_framework.test import APIClient

from classification.enums.classification_enums import ClinicalSignificance, ShareLevel
from classification.tests.models.test_utils import ClassificationTestUtils
from mme.metrics import build_metrics, get_metrics, MME_METRICS_CACHE_KEY
from mme.models import MMEInboundQuery, MMEInboundMatch
from mme.tests.fakes import make_classification

TOKEN = "inbound-secret"


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN},
                   MME_DISCLAIMER="ours", MME_TERMS="our terms")
class MetricsTestCase(TestCase):

    def setUp(self):
        cache.delete(MME_METRICS_CACHE_KEY)
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()

    def tearDown(self):
        cache.delete(MME_METRICS_CACHE_KEY)

    def test_counts_eligible_cases_and_their_genes(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        make_classification(self.lab, self.user, gene_symbol="TP53")

        metrics = build_metrics()
        self.assertEqual(metrics["numberOfCases"], 3)
        self.assertEqual(metrics["numberOfSubmitters"], 1)     # counted as distinct labs
        self.assertEqual(metrics["numberOfGenes"], 3)
        self.assertEqual(metrics["numberOfUniqueGenes"], 2)
        self.assertIn("dateGenerated", metrics)

    def test_ineligible_records_are_excluded(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        make_classification(self.lab, self.user, gene_symbol="TP53",
                            clinical_significance=ClinicalSignificance.BENIGN)
        make_classification(self.lab, self.user, gene_symbol="MSH2",
                            share_level=ShareLevel.ALL_USERS.value)
        make_classification(self.lab, self.user, gene_symbol="MLH1", withdrawn=True)

        self.assertEqual(build_metrics()["numberOfCases"], 1)

    def test_audit_rows_drive_the_traffic_counts(self):
        classification = make_classification(self.lab, self.user, gene_symbol="BRCA1")
        for _ in range(2):
            inbound_query = MMEInboundQuery.objects.create(
                peer_node_id="testpeer", request_json={}, num_results=1)
            MMEInboundMatch.objects.create(
                inbound_query=inbound_query, classification=classification, score=0.5,
                remote_patient_id="remote-1", query_patient_json={})

        metrics = build_metrics()
        self.assertEqual(metrics["numberOfRequestsReceived"], 2)
        self.assertEqual(metrics["numberOfPotentialMatchesSent"], 2)
        self.assertEqual(metrics["numberOfUniqueGenesMatched"], 1)

    def test_get_metrics_serves_cache_until_refreshed(self):
        make_classification(self.lab, self.user, gene_symbol="BRCA1")
        self.assertEqual(get_metrics()["numberOfCases"], 1)

        make_classification(self.lab, self.user, gene_symbol="TP53")
        self.assertEqual(get_metrics()["numberOfCases"], 1)              # still cached
        self.assertEqual(get_metrics(refresh=True)["numberOfCases"], 2)  # recomputed
        self.assertEqual(get_metrics()["numberOfCases"], 2)              # and re-cached


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS={"testpeer": TOKEN},
                   MME_DISCLAIMER="ours", MME_TERMS="our terms")
class MetricsEndpointTestCase(TestCase):

    def setUp(self):
        cache.delete(MME_METRICS_CACHE_KEY)
        self.api_client = APIClient()

    def tearDown(self):
        cache.delete(MME_METRICS_CACHE_KEY)

    def test_api_requires_a_peer_token(self):
        response = self.api_client.get(reverse("mme_api_metrics"))
        self.assertIn(response.status_code, (401, 403))

    def test_api_returns_metrics_and_our_disclaimer(self):
        response = self.api_client.get(reverse("mme_api_metrics"), HTTP_X_AUTH_TOKEN=TOKEN)
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("numberOfCases", data["metrics"])
        self.assertEqual(data["disclaimer"], "ours")
        self.assertEqual(data["terms"], "our terms")

