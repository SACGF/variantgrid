"""
Tests for the node grid reporting a row count it already knows (issue #1257).

The node load pipeline counts each node version, so the grid hands that number to the paginator
rather than making it run a COUNT(*) over the annotated grid queryset.
"""
from django.db import connection
from django.test.utils import CaptureQueriesContext

from analysis.grids import VariantGrid
from analysis.models.nodes.analysis_node import NodeCount
from analysis.tests.test_grid_export import GridExportTestCase
from snpdb.models.models_enums import BuiltInFilters


class NodeGridKnownCountTest(GridExportTestCase):
    FILTERED_POSITION = 2500  # 4 of the 7 test variants are below this

    def _get_data(self, grid, request) -> tuple[dict, list[str]]:
        """ (grid data, the COUNT(*) queries it ran) """
        with CaptureQueriesContext(connection) as queries:
            data = grid.get_data(request)
        count_queries = [q["sql"] for q in queries.captured_queries if q["sql"].startswith("SELECT COUNT(")]
        return data, count_queries

    def test_stored_node_count_skips_the_count_query(self):
        node = self._sample_node()
        grid = VariantGrid(self.user, node)
        data, count_queries = self._get_data(grid, self._request())
        self.assertEqual(node.count, data["records"])
        self.assertEqual([], count_queries)

    def test_extra_filter_uses_its_node_count(self):
        """ NodeCount is keyed on (node_version, label), so the stored count already matches the filter """
        node = self._sample_node()
        NodeCount.objects.create(node_version=node.node_version, label=BuiltInFilters.CLINVAR, count=42)
        grid = VariantGrid(self.user, node, extra_filters=BuiltInFilters.CLINVAR)
        data, count_queries = self._get_data(grid, self._request())
        self.assertEqual(42, data["records"])
        self.assertEqual([], count_queries)

    def test_column_filter_counts_the_queryset(self):
        """ A column filter narrows the rows the stored count was taken over """
        node = self._sample_node()
        grid = VariantGrid(self.user, node)
        data, count_queries = self._get_data(grid, self._request(position_less_than=self.FILTERED_POSITION))
        self.assertLess(data["records"], node.count)
        self.assertEqual(len(data["rows"]), data["records"])
        self.assertEqual(1, len(count_queries))

    def test_node_without_a_count_falls_back_to_the_queryset(self):
        node = self._sample_node()
        expected_records = node.count
        node.count = None
        node.save()

        grid = VariantGrid(self.user, node)
        data, count_queries = self._get_data(grid, self._request())
        self.assertEqual(expected_records, data["records"])
        self.assertEqual(1, len(count_queries))
