from unittest.mock import MagicMock, patch

import requests
from django.test import TestCase

from classification.models.condition_text_search import condition_text_search


def _mock_response(*, json_data=None, json_exc=None, raise_for_status_exc=None) -> MagicMock:
    response = MagicMock()
    if raise_for_status_exc is not None:
        response.raise_for_status.side_effect = raise_for_status_exc
    else:
        response.raise_for_status.return_value = None
    if json_exc is not None:
        response.json.side_effect = json_exc
    else:
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
    def test_non_json_response_returns_empty(self, mock_get):
        """ Rollbar 7269 - Monarch intermittently returns a non-JSON error page, so .json() raises.
        Treat it as no results rather than letting JSONDecodeError escape. """
        mock_get.return_value = _mock_response(json_exc=ValueError("Expecting value: line 2 column 1 (char 1)"))

        self.assertEqual([], condition_text_search("anything"))

    @patch("classification.models.condition_text_search.requests.get")
    def test_http_error_response_returns_empty(self, mock_get):
        """ A 5xx from Monarch (e.g. 502) must not abort the condition search. """
        mock_get.return_value = _mock_response(raise_for_status_exc=requests.exceptions.HTTPError("502 Server Error"))

        self.assertEqual([], condition_text_search("anything"))

    @patch("classification.models.condition_text_search.requests.get")
    def test_connection_error_returns_empty(self, mock_get):
        """ Network failure reaching Monarch must not abort the condition search. """
        mock_get.side_effect = requests.exceptions.ConnectionError("could not connect")

        self.assertEqual([], condition_text_search("anything"))
