"""
Tests for the analysis node grid's two DataTables endpoints (plan Stage 4 item 1).

node_grid_config answers the table definition (columns, renderers, and the node's per-request state
as postData); node_grid_handler answers the rows in the DataTables envelope. Both stay behind
@cache_page, so what the client sends has to be a stable, whitelisted param set.
"""
import json
from urllib.parse import parse_qsl, urlencode

from django.test.client import Client
from django.urls.base import reverse

from analysis.tests.test_grid_export import GridExportTestCase
from analysis.views.views_grid import _NODE_GRID_ALLOWED_PARAMS, _add_allowed_node_grid_params


class NodeGridEndpointTestCase(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.client = Client()
        self.client.force_login(self.user)

    def _config(self) -> dict:
        url = reverse("node_grid_config", kwargs={"analysis_id": self.analysis.pk,
                                                  "analysis_version": self.analysis.version,
                                                  "node_id": self.node.pk,
                                                  "node_version": self.node.version,
                                                  "extra_filters": "default"})
        response = self.client.get(url)
        self.assertEqual(200, response.status_code)
        return json.loads(response.content)

    def _handler_params(self, **extra) -> dict:
        params = {
            "node_id": str(self.node.pk),
            "version_id": str(self.node.version),
            "start": "0",
            "length": "10",
        }
        params.update(extra)
        return params

    def _handler(self, **extra) -> dict:
        url = reverse("node_grid_handler", kwargs={"analysis_id": self.analysis.pk})
        response = self.client.get(url, self._handler_params(**extra))
        self.assertEqual(200, response.status_code)
        return json.loads(response.content)


class NodeGridConfigTest(NodeGridEndpointTestCase):

    def test_definition_shape(self):
        definition = self._config()
        self.assertTrue(definition["columns"])
        self.assertIn("id", {c["data"] for c in definition["columns"]})
        self.assertEqual("GET", definition["ajaxType"])
        self.assertTrue(definition["cacheStableParams"])

    def test_rows_are_deferred(self):
        """ Whether to fetch rows at all is the page's call - the placeholder on a big node, and the
            grid tab being hidden, both hold the query back """
        self.assertTrue(self._config()["deferLoading"])

    def test_no_filter_grid_button_but_the_editor_still_gets_its_fields(self):
        """ Filtering an analysis is what FilterNode is for, so the grid shows no "Filter grid..."
            button - but its editor mounts a builder off the same definition """
        definition = self._config()
        self.assertFalse(definition["filterBuilderToolbar"])
        self.assertTrue(definition["filterBuilder"]["fields"])
        self.assertTrue(definition["filterBuilder"]["operations"])

    def test_post_data_carries_the_nodes_per_request_state(self):
        """ This is what the page sends as its ajax params on every row request """
        post_data = self._config()["postData"]
        self.assertEqual(self.node.pk, post_data["node_id"])
        self.assertEqual(self.node.version, post_data["version_id"])
        ccc = self.analysis.custom_columns_collection
        self.assertEqual(ccc.pk, post_data["ccc_id"])
        self.assertEqual(ccc.version_id, post_data["ccc_version_id"])
        self.assertIn("zygosity_samples_hash", post_data)  # the node has a sample


class NodeGridHandlerTest(NodeGridEndpointTestCase):

    def test_envelope_shape(self):
        data = self._handler()
        self.assertEqual(self.node.count, data["recordsTotal"])
        self.assertEqual(self.node.count, data["recordsFiltered"])
        self.assertEqual(self.node.count, len(data["data"]))
        # Rows keep their .values() field name keys, so they match the columns' 'data' keys
        self.assertIn("locus__position", data["data"][0])

    def test_paging_translates_to_the_grid_engines(self):
        page_two = self._handler(start="2", length="2")
        self.assertEqual(2, len(page_two["data"]))

    def test_column_filters_narrow_the_rows(self):
        """ The filter builder's rules ride up as _search/filters and the grid engine parses them """
        rules = {"groupOp": "AND", "rules": [{"field": "locus__position", "op": "lt", "data": "1500"}]}
        data = self._handler(_search="true", filters=json.dumps(rules))
        self.assertLess(data["recordsFiltered"], self.node.count)


class NodeGridCacheKeyTest(NodeGridEndpointTestCase):
    """ node_grid_handler is @cache_page(WEEK_SECS) - a param that varies per request silently costs
        every user the week-long cache, so the client sends a minimal, stably ordered set """

    # What DataTableDefinition.buildAjaxData emits for this grid, key sorted (@see cacheStableParams)
    CLIENT_PARAMS = [
        ("ccc_id", "1"),
        ("ccc_version_id", "2"),
        ("extra_filters", ""),
        ("length", "10"),
        ("node_id", "3"),
        ("order[0][column]", "4"),
        ("order[0][dir]", "asc"),
        ("start", "0"),
        ("version_id", "5"),
        ("zygosity_samples_hash", "abc123"),
    ]

    def test_every_param_the_client_sends_is_allowed(self):
        """ A param that isn't whitelisted is dropped from the redirect URL, so the retry would ask
            for something other than what the user is looking at """
        for key, _value in self.CLIENT_PARAMS:
            self.assertIn(key, _NODE_GRID_ALLOWED_PARAMS, f"{key} would be dropped")

    def test_identical_grid_state_gives_an_identical_url(self):
        url = "/node_grid/handler/"
        first = _add_allowed_node_grid_params(url, dict(self.CLIENT_PARAMS))
        second = _add_allowed_node_grid_params(url, dict(self.CLIENT_PARAMS))
        self.assertEqual(first, second)
        self.assertEqual(f"{url}?" + urlencode(self.CLIENT_PARAMS), first)

    def test_draw_never_reaches_the_cache_key(self):
        """ 'draw' counts up per request - the client strips it and restores it on the response """
        self.assertNotIn("draw", _NODE_GRID_ALLOWED_PARAMS)
        url = _add_allowed_node_grid_params("/node_grid/handler/",
                                            dict(self.CLIENT_PARAMS + [("draw", "7")]))
        self.assertNotIn("draw", dict(parse_qsl(url.split("?")[1])))
