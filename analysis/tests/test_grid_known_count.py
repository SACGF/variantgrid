"""
Tests for the node grid reporting a row count it already knows (issue #1257).

The node load pipeline counts each node version, so the grid hands that number to the view rather
than making it run a COUNT(*) over the annotated grid queryset.
"""
from django.db import connection
from django.test.utils import CaptureQueriesContext

from analysis.grids import VariantGrid
from analysis.models.nodes.analysis_node import NodeCount
from analysis.tests.test_grid_export import GridExportTestCase
from snpdb.models.models_enums import BuiltInFilters
from snpdb.views.datatable_view import datatable_data


class NodeGridKnownCountTest(GridExportTestCase):
    FILTERED_POSITION = 2500  # 4 of the 7 test variants are below this

    def _get_data(self, grid) -> tuple[dict, list[str]]:
        """ (grid data, the COUNT(*) queries it ran) """
        with CaptureQueriesContext(connection) as queries:
            data = datatable_data(grid)
        count_queries = [q["sql"] for q in queries.captured_queries if q["sql"].startswith("SELECT COUNT(")]
        return data, count_queries

    def _grid(self, request=None, **kwargs) -> VariantGrid:
        if request is None:
            request = self._request()
        return VariantGrid(request, self._sample_node(), **kwargs)

    def test_stored_node_count_skips_the_count_query(self):
        grid = self._grid()
        data, count_queries = self._get_data(grid)
        self.assertEqual(grid.node.count, data["recordsTotal"])
        self.assertEqual(grid.node.count, data["recordsFiltered"])
        self.assertEqual([], count_queries)

    def test_extra_filter_uses_its_node_count(self):
        """ NodeCount is keyed on (node_version, label), so the stored count already matches the filter """
        node = self._sample_node()
        NodeCount.objects.create(node_version=node.node_version, label=BuiltInFilters.CLINVAR, count=42)
        grid = VariantGrid(self._request(), node, extra_filters=BuiltInFilters.CLINVAR)
        data, count_queries = self._get_data(grid)
        self.assertEqual(42, data["recordsTotal"])
        self.assertEqual([], count_queries)

    def test_column_filter_counts_the_queryset(self):
        """ A column filter narrows the rows the stored count was taken over """
        grid = self._grid(request=self._request(position_less_than=self.FILTERED_POSITION))
        data, count_queries = self._get_data(grid)
        self.assertLess(data["recordsFiltered"], grid.node.count)
        self.assertEqual(len(data["data"]), data["recordsFiltered"])
        # Only the filtered count - the whole-table total costs more than the page it's counting
        self.assertEqual(1, len(count_queries))

    def test_node_without_a_count_falls_back_to_the_queryset(self):
        node = self._sample_node()
        expected_records = node.count
        node.count = None
        node.save()

        grid = VariantGrid(self._request(), node)
        data, count_queries = self._get_data(grid)
        self.assertEqual(expected_records, data["recordsTotal"])
        self.assertEqual(1, len(count_queries))
