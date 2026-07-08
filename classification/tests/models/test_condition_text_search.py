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
