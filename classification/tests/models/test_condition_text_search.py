from unittest.mock import MagicMock, patch

from django.test import TestCase

from classification.models.condition_text_search import condition_text_search


def _mock_response(*, json_data=None) -> MagicMock:
    response = MagicMock()
    response.json.return_value = json_data
    return response


class ConditionTextSearchTestCase(TestCase):

    @patch("classification.models.condition_text_search.requests.get")
    def test_unsupported_ontology_prefix_is_skipped(self, mock_get):
        """ #1603 - a result whose prefix isn't an OntologyService (e.g. MPATH) must not abort the
        whole search; supported results around it are still returned. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "MPATH:0000001"},
            {"id": "MONDO:0000001"},
        ]})

        terms = condition_text_search("anything")

        self.assertEqual(["MONDO:0000001"], [term.id for term in terms])

    @patch("classification.models.condition_text_search.requests.get")
    def test_aliased_prefix_is_kept(self, mock_get):
        """ A supported ontology under an alias spelling (Orphanet -> ORPHA) must not be skipped. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "Orphanet:12345"},
        ]})

        terms = condition_text_search("anything")

        self.assertEqual(["ORPHA:12345"], [term.id for term in terms])

    @patch("classification.models.condition_text_search.requests.get")
    def test_malformed_id_is_not_silently_skipped(self, mock_get):
        """ A malformed id is a genuine problem, not an unsupported ontology - it should surface
        (get_or_stub raises) rather than being swallowed like an MPATH skip. """
        mock_get.return_value = _mock_response(json_data={"items": [
            {"id": "not-a-valid-id"},
        ]})

        with self.assertRaises(ValueError):
            condition_text_search("anything")
