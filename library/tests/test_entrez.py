from unittest import mock

from django.test import SimpleTestCase

from library.entrez import EUTILS_TIMEOUT_SECS, eutils_request


class EutilsRequestTest(SimpleTestCase):
    @mock.patch("library.entrez._session")
    def test_sends_timeout_and_comma_joined_ids(self, session):
        eutils_request("efetch", db="pubmed", id=["1", "2"])
        (url,), kwargs = session.post.call_args
        self.assertEqual(url, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi")
        self.assertEqual(kwargs["timeout"], EUTILS_TIMEOUT_SECS)
        self.assertEqual(kwargs["data"]["id"], "1,2")
        session.post.return_value.raise_for_status.assert_called_once()
