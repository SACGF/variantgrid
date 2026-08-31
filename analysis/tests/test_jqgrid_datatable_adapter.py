"""
Tests for the jqGrid -> DataTables adapter (plan Stage 2).

The grid classes are the engine and keep their own tests - these cover the translation layer:
the definition JSON, DataTables params in / jqGrid params out (including the packed genotype sort
index round trip), the data envelope, and the draw handling that keeps the data URL cacheable.
"""
import json

from django.db import connection
from django.test import RequestFactory
from django.test.utils import CaptureQueriesContext
from django.urls.base import reverse

from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils.jqgrid_datatable_adapter import (
    DEFAULT_COLUMN_WIDTH,
    JqGridDatatableView,
    datatable_columns_from_colmodels,
    datatable_filter_fields_from_colmodels,
    datatable_filter_operations,
    datatable_order_from_config,
)
from library.jqgrid.jqgrid import FILTER_OPERATIONS
from library.django_utils.jqgrid_view import JQGridViewOp
from snpdb.models import UserGridConfig
from variantopedia.grids import AllVariantsGrid


class AdapterColumnsTest(GridExportTestCase):
    """ colmodels -> DataTables columns """

    def test_column_shape(self):
        colmodels = [
            {"name": "id", "label": "Variant", "formatter": "representativeVariant", "width": 280,
             "headerTitle": "The variant"},
            {"name": "tags_global", "label": "Tags", "classes": "no-word-wrap",
             "formatter": "tagsGlobalFormatter", "sortable": False},
            {"name": "secret", "label": "Secret", "hidden": True},
        ]
        columns = datatable_columns_from_colmodels(colmodels)

        self.assertEqual(["id", "tags_global", "secret"], [c["data"] for c in columns])
        self.assertEqual("VariantGridFormat.representativeVariant", columns[0]["render"])
        self.assertEqual("The variant", columns[0]["headerTitle"])
        self.assertTrue(columns[0]["orderable"])  # jqGrid colmodels are sortable unless they say otherwise

        self.assertFalse(columns[1]["orderable"])
        self.assertIn("no-word-wrap", columns[1]["className"])

        # Hidden columns stay in the list so a client column index maps back to its colmodel
        self.assertFalse(columns[2]["visible"])
        self.assertIsNone(columns[2]["render"])

    def test_composite_column_extras(self):
        """ A composite cell's renderer kwargs and its header's alternative sort keys """
        colmodels = [
            {"name": "sample_1_samples_zygosity", "formatter": "sampleZygosityFormatter",
             "formatter_kwargs": {"samplePrefix": "sample_1_"},
             "sort_menu": [{"label": "Read depth", "column": "sample_1_samples_read_depth"}]},
            {"name": "sample_1_samples_read_depth", "hidden": True},
        ]
        columns = datatable_columns_from_colmodels(colmodels)
        self.assertEqual({"samplePrefix": "sample_1_"}, columns[0]["renderKwargs"])
        self.assertEqual([{"label": "Read depth", "column": "sample_1_samples_read_depth"}],
                         columns[0]["sortMenu"])
        self.assertNotIn("sortMenu", columns[1])

    def test_every_column_carries_a_width(self):
        """ The client lays these tables out table-layout: fixed, which needs a width per column -
            without one a cell holding 40 PubMed links runs to a few thousand pixels """
        columns = datatable_columns_from_colmodels([{"name": "id", "width": 110}, {"name": "pubmed"}])
        self.assertEqual("110px", columns[0]["width"])
        self.assertEqual(f"{DEFAULT_COLUMN_WIDTH}px", columns[1]["width"])

    def test_order_from_sortname(self):
        colmodels = [{"name": "id"}, {"name": "locus__position", "index": "locus__position"}]
        config = {"sortname": "locus__position", "sortorder": "desc"}
        self.assertEqual([[1, "desc"]], datatable_order_from_config(config, colmodels))

    def test_no_order_when_the_sort_column_is_not_sortable(self):
        """ Sorting is disabled above ANALYSIS_GRID_SORT_MAX_ROWS - an initial order would run it anyway """
        colmodels = [{"name": "locus__position", "sortable": False}]
        config = {"sortname": "locus__position", "sortorder": "asc"}
        self.assertIsNone(datatable_order_from_config(config, colmodels))

    def test_no_order_for_the_default_pk_sortname(self):
        self.assertIsNone(datatable_order_from_config({"sortname": "pk"}, [{"name": "id"}]))


class AdapterParamTranslationTest(GridExportTestCase):

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self.user, self.node)
        self.view = JqGridDatatableView()

    def _translate(self, **params):
        request = RequestFactory().get("/grid/", params)
        request.user = self.user
        return self.view.translate_params(request, self.grid)

    def test_start_and_length_become_page_and_rows(self):
        params = self._translate(start=50, length=25)
        self.assertEqual(25, params["rows"])
        self.assertEqual(3, params["page"])

    def test_first_page(self):
        params = self._translate(start=0, length=10)
        self.assertEqual(1, params["page"])

    def test_no_length_means_no_pagination(self):
        """ rows=0 is the grid engine's 'everything' - the export path relies on it """
        params = self._translate(start=0, length=0)
        self.assertEqual(0, params["rows"])

    def test_a_request_without_datatables_params_pages_off_the_grid_config(self):
        """ A bookmarked grid URL has no 'length' - leave it out so get_paginate_by uses rowNum """
        params = self._translate()
        self.assertNotIn("rows", params)
        self.assertNotIn("page", params)

    def test_order_column_index_becomes_sidx(self):
        colmodels = self.grid.get_colmodels()
        index = next(i for i, cm in enumerate(colmodels) if cm["name"] == "locus__position")
        params = self._translate(**{"start": 0, "length": 10,
                                    "order[0][column]": index, "order[0][dir]": "desc"})
        self.assertEqual("locus__position", params["sidx"])
        self.assertEqual("desc", params["sord"])

    def test_packed_genotype_sort_index_round_trips(self):
        """ Genotype columns pack their sort info into the colmodel 'index'
            (e.g. 'cohortgenotype_134:1:samples_zygosity') - VariantGrid._sort_items decodes it """
        colmodels = self.grid.get_colmodels()
        index, cm = next((i, cm) for i, cm in enumerate(colmodels)
                         if ":" in str(cm.get("index", "")))
        params = self._translate(**{"start": 0, "length": 10, "order[0][column]": index})
        self.assertEqual(cm["index"], params["sidx"])
        self.assertEqual(3, len(params["sidx"].split(":")))

        # ...and the grid turns that into a real sort rather than dropping it
        items = self.grid._sort_items(self.grid.get_queryset(None), params["sidx"], "asc")
        self.assertIn("cohort_genotype_sample_sort_alias", str(items.query))

    def test_column_filters_pass_through_untouched(self):
        filters = json.dumps({"groupOp": "AND", "rules": [{"op": "lt", "field": "locus__position", "data": "2500"}]})
        params = self._translate(start=0, length=10, _search="true", filters=filters)
        self.assertEqual("true", params["_search"])
        self.assertEqual(filters, params["filters"])


class AdapterDataEnvelopeTest(GridExportTestCase):

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self.user, self.node)
        self.view = JqGridDatatableView()

    def _get_data(self, **params):
        request = RequestFactory().get("/grid/", params)
        request.user = self.user
        return self.view.get_data(request, self.grid)

    def test_envelope_shape(self):
        data = self._get_data(start=0, length=10)
        self.assertEqual(self.node.count, data["recordsTotal"])
        self.assertEqual(self.node.count, data["recordsFiltered"])
        self.assertEqual(self.node.count, len(data["data"]))
        # Rows keep their .values() field name keys, so they match the columns' 'data' keys
        self.assertIn("id", data["data"][0])
        self.assertIn("locus__position", data["data"][0])

    def test_paging(self):
        data = self._get_data(start=2, length=2)
        self.assertEqual(2, len(data["data"]))
        self.assertEqual(self.node.count, data["recordsTotal"])

    def test_known_count_runs_no_count_query(self):
        """ The stored node count flows through the adapter - the paginator never runs COUNT(*) """
        with CaptureQueriesContext(connection) as queries:
            data = self._get_data(start=0, length=10)
        count_queries = [q["sql"] for q in queries.captured_queries if q["sql"].startswith("SELECT COUNT(")]
        self.assertEqual([], count_queries)
        self.assertEqual(self.node.count, data["recordsTotal"])

    def test_draw_is_only_echoed_when_it_was_sent(self):
        """ The client strips 'draw' so the request URL stays cacheable, and restores it on the
            response. Echoing a 0 here would look stale to DataTables and the draw would be dropped """
        self.assertNotIn("draw", self._get_data(start=0, length=10))
        self.assertEqual(4, self._get_data(start=0, length=10, draw=4)["draw"])


class AdapterDefinitionTest(GridExportTestCase):

    def _definition(self, **params) -> dict:
        request = RequestFactory().get("/grid/", {"dataTableDefinition": 1, **params})
        request.user = self.user
        grid = VariantGrid(self.user, self._sample_node())
        return JqGridDatatableView().json_definition(request, grid)

    def test_definition_shape(self):
        definition = self._definition()
        self.assertEqual("GET", definition["ajaxType"])  # so @cache_page on the data endpoint works
        self.assertTrue(definition["cacheStableParams"])
        self.assertEqual("VariantGrid", definition["gridName"])
        self.assertTrue(definition["columns"])
        column_names = {c["data"] for c in definition["columns"]}
        self.assertIn("id", column_names)

    def test_analysis_node_visibility_is_grid_wide_not_per_column(self):
        extra = self._definition()["extra"]
        self.assertTrue(extra["analysisNode"]["visible"])
        self.assertEqual(self.genome_build.name, extra["genomeBuild"])

    def test_download_url_reverses_onto_the_download_op(self):
        """ The client renders the toolbar CSV button from it - server side streaming, every row """
        url = reverse("all_variants_grid",
                      kwargs={"genome_build_name": self.genome_build.name, "op": JQGridViewOp.HANDLER})
        request = RequestFactory().get(url, {"dataTableDefinition": 1})
        request.user = self.user
        grid = AllVariantsGrid(self.user, self.genome_build.name)
        definition = JqGridDatatableView(csv_download=True).json_definition(
            request, grid, genome_build_name=self.genome_build.name, op=JQGridViewOp.HANDLER)

        self.assertEqual(reverse("all_variants_grid",
                                 kwargs={"genome_build_name": self.genome_build.name,
                                         "op": JQGridViewOp.DOWNLOAD}),
                         definition["downloadUrl"])
        self.assertFalse(definition["downloadCsvButtonEnabled"], "Client side CSV pages through the grid")

    def test_page_length_comes_from_the_users_grid_config(self):
        UserGridConfig.objects.update_or_create(user=self.user, grid_name="VariantGrid",
                                                defaults={"rows": 15})
        definition = self._definition()
        self.assertEqual(15, definition["pageLength"])
        self.assertIn(15, definition["lengthMenu"])


class AdapterFilterBuilderTest(GridExportTestCase):
    """ The filter builder's fields and operations (plan Stage 4 item 6). The rules it emits go up
        as '_search'/'filters' and are parsed by JqGrid.get_q - the same JSON FilterNode persists. """

    def test_operations_match_the_engines(self):
        operations = datatable_filter_operations()
        self.assertEqual([op for op, _label, _takes in FILTER_OPERATIONS],
                         [o["op"] for o in operations])
        takes_data = {o["op"]: o["takesData"] for o in operations}
        self.assertFalse(takes_data["nu"], "'is null' filters on its own")
        self.assertTrue(takes_data["lt"])

    def test_unfilterable_columns_are_left_out(self):
        """ search=False is jqGrid's own 'don't offer this', and a packed genotype column's 'index'
            is a sample position rather than a field path - neither is a Django lookup """
        fields = datatable_filter_fields_from_colmodels([
            {"name": "locus__position", "label": "Position", "sorttype": "int"},
            {"name": "some_fk", "label": "FK", "search": False},
            {"name": "sample_1_samples_zygosity", "label": "Zyg",
             "index": "cohortgenotype_134:1:samples_zygosity"},
        ])
        self.assertEqual(["locus__position"], [f["field"] for f in fields])
        self.assertEqual("int", fields[0]["type"])

    def test_select_columns_carry_their_choices(self):
        fields = datatable_filter_fields_from_colmodels([
            {"name": "variantannotation__impact", "label": "Impact", "stype": "select",
             "searchoptions": {"value": {"1": "MODIFIER", "4": "HIGH"}}},
        ])
        self.assertEqual("select", fields[0]["type"])
        self.assertEqual({"1": "MODIFIER", "4": "HIGH"}, fields[0]["choices"])

    def test_definition_offers_the_grids_own_columns(self):
        request = RequestFactory().get("/grid/", {"dataTableDefinition": 1})
        request.user = self.user
        grid = VariantGrid(self.user, self._sample_node())
        definition = JqGridDatatableView().json_definition(request, grid)

        fields = definition["filterBuilder"]["fields"]
        offered = {f["field"] for f in fields}
        self.assertIn("locus__position", offered)
        self.assertFalse([f for f in fields if ":" in f["field"]], "Packed genotype columns can't be filtered")

    def test_an_emitted_rule_filters_the_grid(self):
        """ The builder's JSON round trips through the engine that FakeFilterGrid uses """
        node = self._sample_node()
        grid = VariantGrid(self.user, node)
        rules = {"groupOp": "AND",
                 "rules": [{"field": "locus__position", "op": "lt", "data": "1500"}]}
        request = RequestFactory().get("/grid/", {"start": 0, "length": 10,
                                                  "_search": "true", "filters": json.dumps(rules)})
        request.user = self.user

        data = JqGridDatatableView().get_data(request, grid)
        self.assertLess(data["recordsFiltered"], node.count)
        self.assertTrue(all(row["locus__position"] < 1500 for row in data["data"]))
