"""Ensembl API failures (an HTML error page, a connection error) are reported, not raised as a 500 (#869),
and an unversioned Ensembl accession resolves to the latest version."""
from unittest import mock

from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase, override_settings
from django.urls import reverse
from requests import ConnectionError as RequestsConnectionError
from requests import HTTPError, Response

from annotation.fake_data import get_fake_annotation_version
from genes.fake_data import create_fake_transcript_version
from genes.models import TranscriptVersionSequenceInfo
from genes.transcript_errors import NoTranscript
from genes.transcript_sequence_retrieval import TranscriptSequenceFetcher
from library.django_utils.unittest_utils import URLTestCase, _make_test_client
from snpdb.models.models_genome import GenomeBuild

ENSEMBL_ERROR_PAGE = "<!doctype html><html><title>Error: 500 | EMBL-EBI</title></html>"
GRCH38 = mock.Mock(spec=GenomeBuild, name="GRCh38")
GRCH38.name = "GRCh38"


def _response(status_code: int, text: str) -> Response:
    r = Response()
    r.status_code = status_code
    r._content = text.encode()
    r.url = "https://rest.ensembl.org/sequence/id/ENST00000000001?type=cdna"
    return r


class EnsemblFetchErrorTest(SimpleTestCase):
    def test_html_error_page_raises_http_error(self):
        with mock.patch("genes.transcript_sequence_retrieval.requests.get",
                        return_value=_response(500, ENSEMBL_ERROR_PAGE)):
            with self.assertRaises(HTTPError):
                TranscriptSequenceFetcher().fetch_ensembl("ENST00000000001.1", genome_build=GRCH38)

    def test_ok_but_not_json_raises_no_transcript(self):
        with mock.patch("genes.transcript_sequence_retrieval.requests.get",
                        return_value=_response(200, ENSEMBL_ERROR_PAGE)):
            with self.assertRaises(NoTranscript):
                TranscriptSequenceFetcher().fetch_ensembl("ENST00000000001.1", genome_build=GRCH38)


@override_settings(TRANSCRIPT_SEQUENCE_FETCHER_CLASS="genes.transcript_sequence_retrieval.TranscriptSequenceFetcher")
class EnsemblUnversionedSequenceInfoTest(TestCase):
    def test_unversioned_accession_stores_latest_version(self):
        ensembl_json = '{"id": "ENST00000000001", "version": 9, "seq": "ACGT", "molecule": "dna"}'
        with mock.patch("genes.transcript_sequence_retrieval.requests.get",
                        return_value=_response(200, ensembl_json)):
            tvi = TranscriptVersionSequenceInfo.get("ENST00000000001")
        self.assertEqual(tvi.version, 9)
        stored_versions = TranscriptVersionSequenceInfo.objects.filter(transcript_id="ENST00000000001")
        self.assertEqual(list(stored_versions.values_list("version", flat=True)), [9])


class ViewTranscriptVersionApiDownTest(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        cls.transcript_version = create_fake_transcript_version(grch37)

    def setUp(self):
        self.client = _make_test_client()
        self.client.force_login(self.user)

    def test_api_connection_error_renders_warning(self):
        url = reverse("view_transcript_version", kwargs={"transcript_id": self.transcript_version.transcript_id,
                                                         "version": self.transcript_version.version})
        with mock.patch.object(TranscriptVersionSequenceInfo, "get",
                               side_effect=RequestsConnectionError("Ensembl unreachable")):
            response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Could not retrieve transcript sequence information")
        self.assertContains(response, "Ensembl unreachable")
