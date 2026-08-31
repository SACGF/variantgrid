"""
Tests for the representative Variant column and the Classifications column (variantgrid_private#686).

The renderers are client side, so what's tested here is what the server has to get right for them:
the hidden row fields riding along, the colmodels, the genomic sort and the expanded row.
"""
from django.test.client import Client
from django.urls.base import reverse
from django.utils.html import escape

from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from library.django_utils.jqgrid_datatable_adapter import datatable_definition
from snpdb.grid_columns.custom_columns import (
    CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS,
    CLASSIFICATIONS_COLUMN_ROW_FIELDS,
    COMPOSITE_COLUMN_ROW_FIELDS,
    VARIANT_COLUMN_ROW_FIELDS,
)
from snpdb.models import CustomColumnsCollection


class RepresentativeVariantColumnTest(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self.user, self.node)

    def _colmodels_by_name(self) -> dict:
        return {cm["name"]: cm for cm in self.grid.get_colmodels()}

    def test_variant_column_sorts_in_genome_build_order(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "asc")
        coords = [(v.locus.contig.name, v.locus.position) for v in qs]
        self.assertEqual(coords, [("1", 1000), ("1", 2000), ("2", 5000), ("2", 9000),
                                  ("3", 3000), ("10", 500), ("10", 1000)])

    def test_variant_column_sort_descending(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "desc")
        self.assertEqual([(v.locus.contig.name, v.locus.position) for v in qs][:2],
                         [("10", 1000), ("10", 500)])

    def test_renderer_fields_ride_along_hidden(self):
        colmodels = self._colmodels_by_name()
        self.assertEqual(colmodels["id"]["formatter"], "representativeVariant")
        for field in (VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS
                      + COMPOSITE_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS):
            self.assertIn(field, colmodels)
        # Not in the default collection, so hidden, and labelled from the catalogue for the CSV header
        self.assertTrue(colmodels["locus__ref__seq"]["hidden"])
        self.assertEqual(colmodels["locus__ref__seq"]["label"], "Reference")
        self.assertEqual(colmodels["alt__seq"]["label"], "Alt")

    def test_rows_carry_renderer_fields(self):
        request = FakeRequest(user=self.user)
        request.GET = {"rows": "1", "page": "1"}
        row = self.grid.get_data(request)["rows"][0]
        for field in VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS + COMPOSITE_COLUMN_ROW_FIELDS:
            self.assertIn(field, row)
        self.assertNotIn("classifications", row)  # renderer-only column, like tags

    def test_classifications_column_is_not_sortable(self):
        cm = self._colmodels_by_name()["classifications"]
        self.assertEqual(cm["formatter"], "classificationsFormatter")
        self.assertFalse(cm["sortable"])

    def test_composite_columns_keep_their_server_side_formatting(self):
        """ gnomAD popmax AF is formatted server side (unit -> percent) so the CSV matches the grid -
            the client renderer only adds the population beside it """
        colmodels = self._colmodels_by_name()
        self.assertEqual(colmodels["variantannotation__consequence"]["formatter"], "impactConsequenceFormatter")
        popmax_af = colmodels["variantannotation__gnomad_popmax_af"]
        self.assertEqual(popmax_af["formatter"], "gnomadPopmaxFormatter")
        self.assertIn("server_side_formatter", popmax_af)

    def test_definition_declares_row_expansion(self):
        definition = datatable_definition(self.grid)
        annotation_version_id = self.node.analysis.annotation_version_id
        self.assertEqual(definition["expandClientRenderer"],
                         f"variantGridRowDetail.bind(null, {annotation_version_id})")
        self.assertFalse(definition["expandPrefetch"])

    def test_system_default_collection(self):
        columns = list(CustomColumnsCollection.get_system_default().customcolumn_set
                       .order_by("sort_order").values_list("column_id", flat=True))
        self.assertEqual(columns[:4], ["variant", "classifications", "tags", "tags_global"])
        # Coordinates now live in the Variant cell; impact and popmax inside their composite cells
        for removed in ["chrom", "position", "ref", "alt", "svlen", "hgvs_g", "impact", "gnomad_popmax"]:
            self.assertNotIn(removed, columns)


class VariantGridRowDetailViewTest(GridExportTestCase):
    def test_row_detail_renders(self):
        client = Client()
        client.force_login(self.user)
        variant = self.variants[0]
        url = reverse("variant_grid_row_detail",
                      kwargs={"variant_id": variant.pk,
                              "annotation_version_id": self.analysis.annotation_version_id})
        response = client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"v {variant.pk}")
        self.assertContains(response, escape(variant.full_string))  # HGVS/coordinate '>' is escaped
