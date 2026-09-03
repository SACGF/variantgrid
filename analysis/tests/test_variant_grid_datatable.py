"""
Tests for the variant grids as DatatableConfigs (issue #1815).

The grid classes keep their own tests for what they select and count - these cover what they hand
the DataTables client: the definition JSON, the filter builder's fields and operations, the rows
envelope, and the server side CSV download.
"""
import json

from django.test import RequestFactory
from django.urls.base import resolve, reverse

from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from library.django_utils.filter_rules import FILTER_OPERATIONS
from snpdb.models import UserGridConfig
from snpdb.views.datatable_view import (
    DATATABLE_CSV_PARAM,
    DatabaseTableView,
    datatable_definition,
    datatable_response,
)
from variantopedia.grids import AllVariantsGrid


class VariantGridDefinitionTest(GridExportTestCase):

    def _definition(self) -> dict:
        grid = VariantGrid(FakeRequest(user=self.user), self._sample_node())
        return datatable_definition(grid)

    def test_definition_shape(self):
        definition = self._definition()
        self.assertEqual("GET", definition["ajaxType"])  # so @cache_page on the data endpoint works
        self.assertTrue(definition["cacheStableParams"])
        self.assertEqual("VariantGrid", definition["gridName"])
        self.assertIn("id", {c["data"] for c in definition["columns"]})

    def test_every_column_carries_a_width(self):
        """ The client lays these tables out table-layout: fixed, which needs a width per column -
            without one a cell holding 40 PubMed links runs to a few thousand pixels """
        columns = self._definition()["columns"]
        self.assertTrue(columns)
        for column in columns:
            self.assertIn("width", column, column["data"])

    def test_analysis_node_visibility_is_grid_wide_not_per_column(self):
        extra = self._definition()["extra"]
        self.assertTrue(extra["analysisNode"]["visible"])
        self.assertEqual(self.genome_build.name, extra["genomeBuild"])

    def test_multi_value_columns_carry_their_separator(self):
        """ The cell splits on it so a multi-value column reads as a list - @see VEPColumnDef.separator """
        columns = {c["data"]: c for c in self._definition()["columns"]}
        consequence = columns["consequence_impact"]
        members = {m["path"]: m for m in consequence["renderKwargs"]["members"]}
        self.assertEqual("&", members["variantannotation__consequence"]["separator"])
        self.assertNotIn("separator", members["variantannotation__impact"])

    def test_rows_are_deferred(self):
        """ Whether to fetch rows at all is the page's call - the placeholder on a big node, and the
            grid tab being hidden, both hold the query back """
        self.assertTrue(self._definition()["deferLoading"])

    def test_post_data_carries_the_nodes_per_request_state(self):
        post_data = self._definition()["postData"]
        ccc = self.analysis.custom_columns_collection
        self.assertEqual(ccc.pk, post_data["ccc_id"])
        self.assertEqual(ccc.version_id, post_data["ccc_version_id"])
        self.assertIn("zygosity_samples_hash", post_data)  # the node has a sample

    def test_page_length_comes_from_the_users_grid_config(self):
        UserGridConfig.objects.update_or_create(user=self.user, grid_name="VariantGrid",
                                                defaults={"rows": 15})
        definition = self._definition()
        self.assertEqual(15, definition["pageLength"])
        self.assertIn(15, definition["lengthMenu"])

    def test_download_url_streams_the_server_side_csv(self):
        """ The client renders the toolbar CSV button from it - server side streaming, every row """
        url = reverse("all_variants_grid", kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, {"dataTableDefinition": 1})
        request.resolver_match = resolve(url)
        request.user = self.user

        view = DatabaseTableView(column_class=AllVariantsGrid)
        view.request = request
        view.config = view.config_for_request(request)
        definition = view.json_definition()

        self.assertEqual(f"{url}?{DATATABLE_CSV_PARAM}=1", definition["downloadUrl"])
        self.assertFalse(definition["downloadCsvButtonEnabled"], "Client side CSV pages through the grid")


class VariantGridDataTest(GridExportTestCase):

    def _get_data(self, **params) -> dict:
        request = FakeRequest(user=self.user)
        request.GET = params
        grid = VariantGrid(request, self._sample_node())
        return datatable_response(grid, draw=params.get("draw"))

    def test_envelope_shape(self):
        node = self._sample_node()
        data = self._get_data(length="10")
        self.assertEqual(node.count, data["recordsTotal"])
        self.assertEqual(node.count, data["recordsFiltered"])
        self.assertEqual(node.count, len(data["data"]))
        # Rows are keyed by column name, which is what the columns' 'data' keys are
        self.assertIn("id", data["data"][0])
        self.assertIn("locus__position", data["data"][0])

    def test_paging(self):
        node = self._sample_node()
        data = self._get_data(start="2", length="2")
        self.assertEqual(2, len(data["data"]))
        self.assertEqual(node.count, data["recordsTotal"])

    def test_draw_is_only_echoed_when_it_was_sent(self):
        """ The client strips 'draw' so the request URL stays cacheable, and restores it on the
            response. Echoing a 0 here would look stale to DataTables and the draw would be dropped """
        self.assertNotIn("draw", self._get_data(length="10"))
        self.assertEqual(4, self._get_data(length="10", draw="4")["draw"])


class VariantGridFilterBuilderTest(GridExportTestCase):
    """ The filter builder's fields and operations. The rules it emits go up as 'filters' and are
        parsed by filter_rules.rules_to_q - the same JSON FilterNode persists. """

    def _grid(self, **params) -> VariantGrid:
        request = FakeRequest(user=self.user)
        request.GET = params
        return VariantGrid(request, self._sample_node())

    def test_operations_match_the_rule_vocabulary(self):
        operations = datatable_definition(self._grid())["filterBuilder"]["operations"]
        self.assertEqual([op for op, _label, _takes in FILTER_OPERATIONS],
                         [o["op"] for o in operations])
        takes_data = {o["op"]: o["takesData"] for o in operations}
        self.assertFalse(takes_data["nu"], "'is null' filters on its own")
        self.assertTrue(takes_data["lt"])

    def test_only_columns_that_are_django_lookups_are_offered(self):
        """ 'field' goes straight into a Django lookup, so display only columns (the packed genotype
            cells, the annotation aliases) and ForeignKey endpoints are left out """
        grid = self._grid()
        fields = datatable_definition(grid)["filterBuilder"]["fields"]
        offered = {f["field"] for f in fields}
        self.assertIn("locus__position", offered)
        self.assertEqual("int", next(f for f in fields if f["field"] == "locus__position")["type"])
        self.assertNotIn("tags", offered)
        self.assertFalse([f for f in offered if f.startswith("sample_")],
                         "Packed genotype columns can't be filtered")

    def test_select_columns_carry_their_choices(self):
        fields = {f["field"]: f for f in datatable_definition(self._grid())["filterBuilder"]["fields"]}
        impact = fields["variantannotation__impact"]
        self.assertEqual("select", impact["type"])
        self.assertIn("HIGH", impact["choices"].values())

    def test_an_emitted_rule_filters_the_grid(self):
        node = self._sample_node()
        rules = {"groupOp": "AND",
                 "rules": [{"field": "locus__position", "op": "lt", "data": "1500"}]}
        grid = self._grid(length="10", filters=json.dumps(rules))
        data = datatable_response(grid)
        self.assertLess(data["recordsFiltered"], node.count)
        self.assertTrue(all(row["locus__position"] < 1500 for row in data["data"]))
