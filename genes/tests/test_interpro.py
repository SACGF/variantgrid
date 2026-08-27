from http import HTTPStatus
from unittest.mock import patch

from django.test import TestCase, override_settings
from requests import HTTPError

from genes import interpro
from genes.models import Pfam, PfamDomains, PfamSequence

SEQ_ID = "P38398"
RING_FINGER = {
    "metadata": {"accession": "PF00097", "name": "Zinc finger, C3HC4 type (RING finger)"},
    "proteins": [{
        "accession": SEQ_ID.lower(),
        "protein_length": 1863,
        "entry_protein_locations": [{"fragments": [{"start": 24, "end": 64, "dc-status": "CONTINUOUS"}]}],
    }],
}


def _entry(accession: str, fragments: list[dict]) -> dict:
    return {
        "metadata": {"accession": accession, "name": accession},
        "proteins": [{"accession": SEQ_ID.lower(), "entry_protein_locations": [{"fragments": fragments}]}],
    }


class MockResponse:
    def __init__(self, status_code: int, data=None):
        self.status_code = status_code
        self._data = data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise HTTPError(f"{self.status_code}")

    def json(self):
        return self._data


def mock_results(*entries) -> MockResponse:
    return MockResponse(HTTPStatus.OK, {"count": len(entries), "next": None, "results": list(entries)})


@override_settings(UNIT_TEST=False, PFAM_INTERPRO_LAZY_DOMAINS=True, PFAM_INTERPRO_MAX_WORKERS=2)
class InterProTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        Pfam.objects.create(pk=97, pfam_id="PF00097", description="Zinc finger, C3HC4 type (RING finger)")

    def test_no_pfam_match_is_stamped_with_no_domains(self):
        """ 204 = the protein genuinely has no Pfam entries - has to be remembered or we re-fetch forever """
        pfam_sequence = PfamSequence.objects.create(seq_id=SEQ_ID)
        with patch.object(interpro._session, "get", return_value=MockResponse(HTTPStatus.NO_CONTENT)):
            self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 0)

        pfam_sequence.refresh_from_db()
        self.assertIsNotNone(pfam_sequence.domains_imported)
        self.assertFalse(pfam_sequence.pfamdomains_set.exists())

    def test_already_imported_sequence_is_not_re_fetched(self):
        PfamSequence.objects.create(seq_id=SEQ_ID, domains_imported="2026-01-01T00:00:00Z")
        with patch.object(interpro._session, "get") as mock_get:
            self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 0)
        mock_get.assert_not_called()

    def test_each_fragment_becomes_a_domain(self):
        """ A discontinuous domain draws as separate boxes """
        PfamSequence.objects.create(seq_id=SEQ_ID)
        fragments = [{"start": 24, "end": 64}, {"start": 120, "end": 180}]
        response = mock_results(_entry("PF00097", fragments))
        with patch.object(interpro._session, "get", return_value=response):
            self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 2)

        starts_ends = list(PfamDomains.objects.order_by("start").values_list("start", "end"))
        self.assertEqual(starts_ends, [(24, 64), (120, 180)])

    def test_unknown_pfam_accession_is_skipped(self):
        PfamSequence.objects.create(seq_id=SEQ_ID)
        response = mock_results(RING_FINGER, _entry("PF99999", [{"start": 200, "end": 300}]))
        with patch.object(interpro._session, "get", return_value=response):
            with self.assertLogs(level="WARNING") as logs:
                self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 1)

        self.assertTrue(any("skipped 1 domains" in message for message in logs.output))
        self.assertEqual(PfamDomains.objects.get().pfam_id, 97)

    def test_failed_request_leaves_sequence_unimported(self):
        pfam_sequence = PfamSequence.objects.create(seq_id=SEQ_ID)
        with patch.object(interpro._session, "get", return_value=MockResponse(HTTPStatus.SERVICE_UNAVAILABLE)):
            self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 0)

        pfam_sequence.refresh_from_db()
        self.assertIsNone(pfam_sequence.domains_imported)

    def test_re_fetch_replaces_existing_domains(self):
        pfam_sequence = PfamSequence.objects.create(seq_id=SEQ_ID)
        PfamDomains.objects.create(pfam_sequence=pfam_sequence, pfam_id=97, start=1, end=10)

        response = mock_results(RING_FINGER)
        with patch.object(interpro._session, "get", return_value=response):
            self.assertEqual(interpro.store_domains_for_sequences([SEQ_ID]), 1)

        self.assertEqual(list(PfamDomains.objects.values_list("start", "end")), [(24, 64)])
