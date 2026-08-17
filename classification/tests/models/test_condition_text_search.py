from unittest.mock import MagicMock, patch

from django.test import TestCase
from requests import HTTPError
from requests.exceptions import JSONDecodeError

from classification.models.condition_text_search import condition_text_search


def _mock_response(*, json_data=None) -> MagicMock:
    response = MagicMock()
    response.json.return_value = json_data
    return response


class ConditionTextSearchTestCase(TestCase):

    @patch("classification.models.condition_text_search._session.get")
    def test_unsupported_ontology_prefix_is_skipped(self, mock_get):
        """ #1603 - a result whose prefix isn't an OntologyService (e.g. MPATH) must not abort the
        whole search; supported results around it are still returned. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "MPATH:0000001"},
            {"id": "MONDO:0000001"},
        ]})

        terms = condition_text_search("anything")

        self.assertEqual(["MONDO:0000001"], [term.id for term in terms])

    @patch("classification.models.condition_text_search._session.get")
    def test_aliased_prefix_is_kept(self, mock_get):
        """ A supported ontology under an alias spelling (Orphanet -> ORPHA) must not be skipped. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "Orphanet:12345"},
        ]})

        terms = condition_text_search("anything")

        self.assertEqual(["ORPHA:12345"], [term.id for term in terms])

    @patch("classification.models.condition_text_search._session.get")
    def test_malformed_id_is_skipped_without_aborting(self, mock_get):
        """ A malformed id from Monarch must not abort the whole search - it's skipped like any other
        id we can't turn into a supported ontology term, and valid results around it survive. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "not-a-valid-id"},
            {"id": "MONDO:0000001"},
        ]})

        terms = condition_text_search("anything")

        self.assertEqual(["MONDO:0000001"], [term.id for term in terms])

    @patch("classification.models.condition_text_search._session.get")
    def test_non_json_body_names_monarch(self, mock_get):
        """ #1742 - a 200 with a non-JSON body has to propagate as an error that identifies Monarch
        and what it sent, rather than a bare JSONDecodeError. """
        response = MagicMock()
        response.json.side_effect = JSONDecodeError("Expecting value", "\n<html>502</html>", 1)
        response.headers = {"Content-Type": "text/html"}
        response.text = "\n<html>502 Bad Gateway</html>"
        mock_get.return_value = response

        with self.assertRaises(ValueError) as cm:
            condition_text_search("anything")

        message = str(cm.exception)
        self.assertIn("Monarch", message)
        self.assertIn("text/html", message)
        self.assertIn("502 Bad Gateway", message)

    @patch("classification.models.condition_text_search._session.get")
    def test_http_error_propagates(self, mock_get):
        """ #1742 - an error status propagates as an HTTPError (which names the status and Monarch's
        URL) instead of dying on the JSON parse. """
        response = MagicMock()
        response.raise_for_status.side_effect = HTTPError("503 Server Error for url: monarch")
        mock_get.return_value = response

        with self.assertRaises(HTTPError):
            condition_text_search("anything")
