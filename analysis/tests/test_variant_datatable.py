"""
Tests for the variant grids as DatatableConfigs (issue #1785 section 3).

The grid classes keep their own behavioural tests (counts, sorting limits, export) - these cover what
the DataTables layer makes of them: the columns built from a CustomColumnsCollection, the definition
JSON, the data envelope, the sort annotation that replaced jqGrid's packed 'index' string, and the
filter builder whose rules FilterNode also persists.
"""
import json

from django.test import RequestFactory
from django.urls.base import reverse

from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils.filter_rules import FILTER_OPERATIONS
from snpdb.grid_columns.variant_columns import (
    DEFAULT_COLUMN_WIDTH,
    build_variant_columns,
    column_by_name,
)
from snpdb.models import UserGridConfig, Variant
from snpdb.views.datatable_view import datatable_data, datatable_definition
from variantopedia.grids import AllVariantsGrid


class VariantColumnsTest(GridExportTestCase):
    """ {field: overrides} -> RichColumns -> the definition's columns """

    OVERRIDES = {
        "id": {"label": "Variant", "client_renderer": "VariantGridFormat.detailsLink", "width": 110,
               "header_title": "The variant"},
        "tags_global": {"model_field": False, "queryset_field": False, "label": "Tags",
                        "css_class": "no-word-wrap", "client_renderer": "VariantGridFormat.tagsGlobal",
                        "orderable": False},
        "locus__position": {"visible": False},
    }

    def _columns(self):
        return build_variant_columns(Variant, list(self.OVERRIDES), self.OVERRIDES)

    def test_column_shape(self):
        columns = self._columns()
        self.assertEqual(["id", "tags_global", "locus__position"], [rc.name for rc in columns])

        variant_id = columns[0]
        self.assertEqual("VariantGridFormat.detailsLink", variant_id.client_renderer)
        self.assertEqual("The variant", variant_id.header_title)
        self.assertTrue(variant_id.orderable)  # columns are orderable unless they say otherwise
        self.assertEqual("int", variant_id.data_type)  # taken off the Django field

        tags_global = columns[1]
        self.assertFalse(tags_global.orderable)
        self.assertIn("no-word-wrap", tags_global.css_classes)
        self.assertFalse(tags_global.queryset_field)

        # Hidden columns stay in the list so a client column index maps back to its column
        self.assertFalse(columns[2].visible)

    def test_every_column_carries_a_width(self):
        """ The client lays these tables out table-layout: fixed, which needs a width per column -
            without one a cell holding 40 PubMed links runs to a few thousand pixels """
        columns = self._columns()
        self.assertEqual(110, columns[0].width)
        self.assertEqual(DEFAULT_COLUMN_WIDTH, columns[1].width)

    def test_definition_columns_carry_the_client_renderer_and_width(self):
        grid = VariantGrid(self._request(), self._sample_node())
        definition = datatable_definition(grid)
        by_name = {c["data"]: c for c in definition["columns"]}
        self.assertEqual("VariantGridFormat.detailsLink", by_name["id"]["render"])
        self.assertTrue(all(c["width"].endswith("px") for c in definition["columns"]))


class VariantGridOrderTest(GridExportTestCase):

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self._request(), self.node)

    def _order_params(self, column_name: str, direction="asc") -> dict:
        index = [rc.name for rc in self.grid.enabled_columns].index(column_name)
        return {"order[0][column]": index, "order[0][dir]": direction}

    def _set_analysis_default_sort(self, variant_column: str):
        ccc = self.analysis.custom_columns_collection
        self.analysis.default_sort_by_column = ccc.customcolumn_set.get(column__variant_column=variant_column)
        self.analysis.save()

    def _ordered(self, column_name: str, direction="asc"):
        request = RequestFactory().get("/grid/", self._order_params(column_name, direction))
        request.user = self.user
        grid = VariantGrid(request, self.node)
        return grid.ordering(grid.get_initial_queryset())

    def test_order_by_a_field_column(self):
        query = str(self._ordered("locus__position", "desc").query)
        self.assertIn("locus", query.lower())
        self.assertIn("DESC", query)

    def test_rows_with_no_value_sort_as_the_lowest(self):
        """ Ascending puts them first, descending last - rather than always at the bottom """
        self.assertIn("NULLS FIRST", str(self._ordered("variantannotation__gnomad_af", "asc").query))
        self.assertIn("NULLS LAST", str(self._ordered("variantannotation__gnomad_af", "desc").query))

    def test_packed_genotype_column_sorts_by_its_annotation(self):
        """ A sample's zygosity lives in a slot of its cohort's packed column - there's no field to
            sort on, so the column carries the expression to annotate and order by """
        packed = next(rc for rc in self.grid.enabled_columns if rc.sort_annotation is not None)
        self.assertTrue(packed.name.startswith("sample_"))
        self.assertFalse(packed.filterable, "A packed genotype slot isn't a Django lookup")

        ordered = self._ordered(packed.name)
        self.assertIn(packed.sort_alias, str(ordered.query))

    def test_initial_order_comes_from_the_analysis_default_sort_column(self):
        """ Set on the column, so it's both the grid's initial order and the secondary sort """
        self.assertIsNone(self.grid.default_sort_order_column, "No analysis default sort column set")
        self.assertIsNone(datatable_definition(self.grid).get("order"))

        self._set_analysis_default_sort("locus__position")

        grid = VariantGrid(self._request(), self.node)
        self.assertEqual("locus__position", grid.default_sort_order_column.name)
        column_index = [rc.name for rc in grid.enabled_columns].index("locus__position")
        self.assertEqual([[column_index, "asc"]], datatable_definition(grid)["order"])

    def test_a_sort_on_another_column_keeps_the_default_as_its_secondary(self):
        self._set_analysis_default_sort("locus__position")

        ordered = self._ordered("id", "desc")
        order_by = [str(o.expression) for o in ordered.query.order_by]
        self.assertEqual(["F(id)", "F(locus__position)", "F(pk)"], order_by)


class VariantGridDataEnvelopeTest(GridExportTestCase):

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()

    def _get_data(self, **params):
        request = RequestFactory().get("/grid/", params)
        request.user = self.user
        return datatable_data(VariantGrid(request, self.node))

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

    def test_page_length_is_not_capped(self):
        """ Rows per page is the user's own UserGridConfig setting, so the 100 row cap is off """
        self.assertIsNone(VariantGrid(self._request(), self.node).max_page_length)

    def test_draw_is_only_echoed_when_it_was_sent(self):
        """ The client strips 'draw' so the request URL stays cacheable, and restores it on the
            response. Echoing a 0 here would look stale to DataTables and the draw would be dropped """
        self.assertNotIn("draw", self._get_data(start=0, length=10))
        self.assertEqual(4, self._get_data(start=0, length=10, draw=4)["draw"])


class VariantGridDefinitionTest(GridExportTestCase):

    def _definition(self, **params) -> dict:
        request = RequestFactory().get("/grid/", {"dataTableDefinition": 1, **params})
        request.user = self.user
        return datatable_definition(VariantGrid(request, self._sample_node()))

    def test_definition_shape(self):
        definition = self._definition()
        self.assertEqual("GET", definition["ajaxType"])  # so @cache_page on the data endpoint works
        self.assertTrue(definition["cacheStableParams"])
        self.assertTrue(definition["deferLoading"])  # the page decides when to fetch rows
        self.assertEqual("VariantGrid", definition["gridName"])
        self.assertTrue(definition["columns"])
        self.assertIn("id", {c["data"] for c in definition["columns"]})
        # A node export can take minutes, so it goes through Celery rather than streaming from here
        self.assertNotIn("downloadUrl", definition)

    def test_analysis_node_visibility_is_grid_wide_not_per_column(self):
        extra = self._definition()["extra"]
        self.assertTrue(extra["analysisNode"]["visible"])
        self.assertEqual(self.genome_build.name, extra["genomeBuild"])

    def test_download_url_is_the_grids_own_url(self):
        """ The client renders the toolbar CSV button from it - server side streaming, every row """
        url = reverse("all_variants_grid", kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, {"dataTableDefinition": 1,
                                             "genome_build_name": self.genome_build.name})
        request.user = self.user
        definition = datatable_definition(AllVariantsGrid(request))

        self.assertEqual(f"{url}?dataTableCsv=1", definition["downloadUrl"])
        self.assertFalse(definition["downloadCsvButtonEnabled"], "Client side CSV pages through the grid")

    def test_page_length_comes_from_the_users_grid_config(self):
        UserGridConfig.objects.update_or_create(user=self.user, grid_name="VariantGrid",
                                                defaults={"rows": 15})
        definition = self._definition()
        self.assertEqual(15, definition["pageLength"])
        self.assertIn(15, definition["lengthMenu"])


class VariantGridFilterBuilderTest(GridExportTestCase):
    """ The filter builder's fields and operations. The rules it emits go up as a 'filters' param and
        are parsed by filter_rules_to_q - the same JSON FilterNode persists """

    def _grid(self, **params) -> VariantGrid:
        request = RequestFactory().get("/grid/", params)
        request.user = self.user
        return VariantGrid(request, self._sample_node())

    def test_operations_match_the_engines(self):
        operations = datatable_definition(self._grid())["filterBuilder"]["operations"]
        self.assertEqual([op for op, _label, _takes in FILTER_OPERATIONS],
                         [o["op"] for o in operations])
        takes_data = {o["op"]: o["takesData"] for o in operations}
        self.assertFalse(takes_data["nu"], "'is null' filters on its own")
        self.assertTrue(takes_data["lt"])

    def test_definition_offers_the_grids_own_filterable_columns(self):
        grid = self._grid()
        fields = datatable_definition(grid)["filterBuilder"]["fields"]
        offered = {f["field"] for f in fields}
        self.assertIn("locus__position", offered)
        self.assertEqual("int", next(f["type"] for f in fields if f["field"] == "locus__position"))

        # A packed genotype slot is a sample position rather than a field path, so it can't be filtered
        packed = next(rc for rc in grid.enabled_columns if rc.sort_annotation is not None)
        self.assertNotIn(packed.name, offered)

    def test_select_columns_carry_their_choices(self):
        """ A column whose Django field has choices gets a dropdown rather than a text input """
        grid = self._grid()
        impact = column_by_name(grid.enabled_columns, "variantannotation__impact")
        if impact is None:
            self.skipTest("No choices column in this deployment's default columns")
        self.assertTrue(impact.filter_choices)

    def test_an_emitted_rule_filters_the_grid(self):
        rules = {"groupOp": "AND", "rules": [{"field": "locus__position", "op": "lt", "data": "1500"}]}
        grid = self._grid(start=0, length=10, filters=json.dumps(rules))
        data = datatable_data(grid)
        self.assertLess(data["recordsFiltered"], grid.node.count)
        self.assertTrue(all(row["locus__position"] < 1500 for row in data["data"]))
