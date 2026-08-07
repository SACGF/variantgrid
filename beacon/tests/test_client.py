from types import SimpleNamespace
from unittest.mock import patch, MagicMock

import requests
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from beacon.client import query_node, query_external_beacons_for_variant
from beacon.models import BeaconQueryCache
from beacon.query_targets import eligible_queries, evaluate_queries
from beacon.variant_mapping import variant_to_beacon_query_params, parse_beacon_response
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

NODES = {
    "nodeA": {"base_url": "https://a.beacon", "api_version": "v2.0.0", "token": None, "type": "snv"},
    "nodeB": {"base_url": "https://b.beacon", "api_version": "v2.0.0", "token": "secret", "type": "snv"},
}


def _mock_response(payload):
    response = MagicMock()
    response.raise_for_status.return_value = None
    response.json.return_value = payload
    return response


@override_settings(BEACON_OUTBOUND_ENABLED=True, BEACON_QUERY_NODES=NODES,
                   BEACON_QUERY_TIMEOUT=5, BEACON_QUERY_CACHE_DAYS=7)
class BeaconClientTestCase(TestCase):

    def setUp(self):
        self.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant = slowly_create_test_variant("3", 1000, "A", "T", self.grch37)

    def test_param_building_and_response_parsing(self):
        params = variant_to_beacon_query_params(self.variant, self.grch37)
        self.assertEqual(params["referenceName"], "3")
        self.assertEqual(params["start"], 999)  # 1-based 1000 -> 0-based
        self.assertEqual(params["referenceBases"], "A")
        self.assertEqual(params["alternateBases"], "T")

        exists, count = parse_beacon_response({"responseSummary": {"exists": True, "numTotalResults": 4}})
        self.assertTrue(exists)
        self.assertEqual(count, 4)

    def test_parse_response_distinguishes_absent_from_unanswered(self):
        """ `exists: null` / a missing summary means "cannot answer", not "absent" - collapsing
            either to False would render a definitive "not found" the node never claimed. """
        exists, _ = parse_beacon_response({"responseSummary": {"exists": False, "numTotalResults": 0}})
        self.assertIs(exists, False)

        for unanswered in ({"responseSummary": {"exists": None, "numTotalResults": 0}},
                           {"responseSummary": None},
                           {"responseSummary": {}},
                           {}):
            exists, _ = parse_beacon_response(unanswered)
            self.assertIsNone(exists, f"expected None for {unanswered}")

    def test_query_node_success(self):
        payload = {"responseSummary": {"exists": True, "numTotalResults": 2}}
        with patch("beacon.client.requests.get", return_value=_mock_response(payload)) as mock_get:
            result = query_node("nodeB", {"referenceName": "3"})
        self.assertTrue(result["exists"])
        self.assertEqual(result["count"], 2)
        self.assertIsNone(result["error"])
        # nodeB carries a token -> Bearer header sent
        _args, kwargs = mock_get.call_args
        self.assertEqual(kwargs["headers"]["Authorization"], "Bearer secret")

    def test_query_node_error_is_isolated(self):
        response = MagicMock()
        response.raise_for_status.side_effect = requests.Timeout("timed out")
        with patch("beacon.client.requests.get", return_value=response):
            result = query_node("nodeA", {"referenceName": "3"})
        self.assertIsNone(result["exists"])  # errored -> unknown, never a negative
        self.assertIn("timed out", result["error"])

    def test_fanout_one_node_down_still_returns_others(self):
        def fake_get(url, **kwargs):
            if url.startswith("https://a.beacon"):
                raise requests.ConnectionError("node A down")
            return _mock_response({"responseSummary": {"exists": True, "numTotalResults": 1}})

        with patch("beacon.client.requests.get", side_effect=fake_get):
            results = {r["node_id"]: r for r in query_external_beacons_for_variant(self.variant, self.grch37)}
        self.assertIsNotNone(results["nodeA"]["error"])
        self.assertTrue(results["nodeB"]["exists"])


@override_settings(BEACON_OUTBOUND_ENABLED=True, BEACON_QUERY_NODES=NODES,
                   BEACON_QUERY_TIMEOUT=5, BEACON_QUERY_CACHE_DAYS=7)
class ExternalBeaconsViewTestCase(TestCase):

    def setUp(self):
        self.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant = slowly_create_test_variant("3", 1000, "A", "T", self.grch37)
        self.user = User.objects.create_user("beacon_view_user")

    def _url(self):
        return reverse("external_beacons", kwargs={"variant_id": self.variant.pk,
                                                    "genome_build_name": "GRCh37"})

    def test_renders_and_caches_then_serves_from_cache(self):
        payload = {"responseSummary": {"exists": True, "numTotalResults": 3}}
        self.client.force_login(self.user)
        with patch("beacon.client.requests.get", return_value=_mock_response(payload)) as mock_get:
            response = self.client.get(self._url())
            self.assertEqual(response.status_code, 200)
            first_call_count = mock_get.call_count
            self.assertEqual(BeaconQueryCache.objects.filter(variant=self.variant).count(), 2)

            # Second view: fresh cache -> no further network calls.
            self.client.get(self._url())
            self.assertEqual(mock_get.call_count, first_call_count)

    @override_settings(BEACON_OUTBOUND_ENABLED=False)
    def test_disabled_returns_empty(self):
        self.client.force_login(self.user)
        response = self.client.get(self._url())
        self.assertEqual(response.status_code, 200)
        self.assertNotIn(b"External Beacons", response.content)

    @override_settings(BEACON_QUERY_NODES={
        "progenetix": {"base_url": "https://progenetix.org/beacon", "type": "cnv",
                       "assemblies": ["GRCh38"]}})
    def test_skipped_node_is_explained(self):
        # setUp's variant is a GRCh37 SNV, so the GRCh38-only CNV node is skipped, not queried.
        self.client.force_login(self.user)
        response = self.client.get(self._url())
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"External Beacons", response.content)
        self.assertIn(b"progenetix", response.content)
        self.assertIn(b"GRCh38", response.content)  # explains it's a GRCh38-only Beacon


class BeaconQueryGatingTestCase(TestCase):
    """ eligible_queries() routes each variant only to servers whose domain it matches. """

    def setUp(self):
        self.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        self.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.both_nodes = {"seq": {"type": "snv"}, "cnv": {"type": "cnv"}}

    @staticmethod
    def _variant(**coordinate):
        return SimpleNamespace(coordinate=SimpleNamespace(**coordinate))

    def test_snv_routes_to_snv_node_only(self):
        variant = self._variant(is_symbolic=False, chrom="3", position=1000,
                                ref="A", alt="T", svlen=None)
        eligible = eligible_queries(variant, self.grch38, self.both_nodes)
        self.assertEqual(list(eligible), ["seq"])
        self.assertEqual(eligible["seq"]["start"], 999)  # 1-based -> 0-based
        self.assertEqual(eligible["seq"]["alternateBases"], "T")

    def test_cnv_deletion_routes_to_cnv_node_as_range_query(self):
        variant = self._variant(is_symbolic=True, chrom="9", position=21967752,
                                ref="N", alt="<DEL>", svlen=-27549, end=21995301)
        eligible = eligible_queries(variant, self.grch38, self.both_nodes)
        self.assertEqual(list(eligible), ["cnv"])
        params = eligible["cnv"]
        self.assertEqual(params["start"], 21967751)  # 0-based
        self.assertEqual(params["end"], 21995301)    # position + abs(svlen)
        self.assertEqual(params["variantType"], "EFO:0030067")  # <DEL> -> copy number loss
        self.assertNotIn("referenceBases", params)   # a range query, not an exact-allele one

    def test_cnv_duplication_maps_to_gain(self):
        variant = self._variant(is_symbolic=True, chrom="2", position=15940550,
                                ref="N", alt="<DUP>", svlen=5000, end=15945550)
        eligible = eligible_queries(variant, self.grch38, {"cnv": {"type": "cnv"}})
        self.assertEqual(eligible["cnv"]["variantType"], "EFO:0030070")  # <DUP> -> gain

    def test_assembly_gate_skips_wrong_build(self):
        variant = self._variant(is_symbolic=True, chrom="9", position=21967752,
                                ref="N", alt="<DEL>", svlen=-27549, end=21995301)
        nodes = {"cnv": {"type": "cnv", "assemblies": ["GRCh38"]}}
        self.assertEqual(eligible_queries(variant, self.grch37, nodes), {})

    def test_unsupported_symbolic_type_is_skipped(self):
        variant = self._variant(is_symbolic=True, chrom="1", position=100,
                                ref="N", alt="<INV>", svlen=500, end=600)
        self.assertEqual(eligible_queries(variant, self.grch38, {"cnv": {"type": "cnv"}}), {})

    def test_unknown_node_type_is_skipped(self):
        variant = self._variant(is_symbolic=False, chrom="3", position=1000,
                                ref="A", alt="T", svlen=None)
        self.assertEqual(eligible_queries(variant, self.grch38, {"x": {"type": "mystery"}}), {})

    def test_evaluate_queries_reports_skip_reasons(self):
        # An SNV against a GRCh38-only CNV node: skipped, and evaluate_queries says why.
        variant = self._variant(is_symbolic=False, chrom="3", position=1000,
                                 ref="A", alt="T", svlen=None)
        nodes = {"cnv": {"type": "cnv", "assemblies": ["GRCh38"]}}

        [wrong_build] = evaluate_queries(variant, self.grch37, nodes)
        self.assertFalse(wrong_build.eligible)
        self.assertIn("GRCh38", wrong_build.reason)

        [wrong_domain] = evaluate_queries(variant, self.grch38, nodes)
        self.assertFalse(wrong_domain.eligible)
        self.assertIn("copy-number", wrong_domain.reason)
