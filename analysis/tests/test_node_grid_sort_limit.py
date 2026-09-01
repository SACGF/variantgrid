from django.contrib.auth.models import User
from django.db.models import F
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.grids import VariantGrid
from analysis.models import AllVariantsNode, Analysis
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils import FakeRequest
from snpdb.models import GenomeBuild


@override_settings(ANALYSIS_GRID_SORT_MAX_ROWS=10_000)
class NodeGridSortLimitTest(TestCase):
    """ Large nodes (count >= ANALYSIS_GRID_SORT_MAX_ROWS) disable sorting and fall back to -pk,
        so sorting by a joined/unindexed column can't blow the statement_timeout (issue #1651). """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_node_grid_sort_limit")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    def _grid(self, count, **params):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=count)
        request = FakeRequest(user=self.user)
        request.GET = params
        return VariantGrid(request, node)

    def test_small_node_sorting_enabled(self):
        grid = self._grid(count=500)
        self.assertFalse(grid.sorting_disabled())

    def test_large_node_sorting_disabled(self):
        grid = self._grid(count=10_000)
        self.assertTrue(grid.sorting_disabled())

    def test_unknown_count_sorting_disabled(self):
        grid = self._grid(count=None)
        self.assertTrue(grid.sorting_disabled())

    def _position_column_index(self, grid) -> int:
        return next(i for i, rc in enumerate(grid.enabled_columns) if rc.name == "locus__position")

    def test_large_node_ignores_requested_sort(self):
        grid = self._grid(count=50_000)
        grid.request.GET = {"order[0][column]": str(self._position_column_index(grid)),
                            "order[0][dir]": "asc"}
        sorted_qs = grid.ordering(grid._get_base_queryset())
        self.assertEqual(list(sorted_qs.query.order_by), ["-pk"])

    def test_small_node_keeps_requested_sort(self):
        grid = self._grid(count=500)
        grid.request.GET = {"order[0][column]": str(self._position_column_index(grid)),
                            "order[0][dir]": "asc"}
        sorted_qs = grid.ordering(grid._get_base_queryset())
        # Requested column first, PK tiebreaker last
        self.assertEqual(str(list(sorted_qs.query.order_by)[-1]), str(F("pk").desc()))
        self.assertGreater(len(sorted_qs.query.order_by), 1)

    def test_large_node_columns_not_orderable(self):
        grid = self._grid(count=50_000)
        self.assertTrue(grid.enabled_columns)
        self.assertTrue(all(rc.orderable is False for rc in grid.enabled_columns))

    def test_small_node_columns_orderable_default(self):
        grid = self._grid(count=500)
        # Sorting isn't force-disabled - columns keep their normal orderable setting
        self.assertTrue(any(rc.orderable for rc in grid.enabled_columns))

    def test_large_node_has_no_initial_order(self):
        grid = self._grid(count=50_000)
        self.assertIsNone(grid.initial_order())


@override_settings(ANALYSIS_GRID_SORT_MAX_ROWS=10_000)
class NodeGridSortLimitBannerTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_node_grid_sort_banner")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    def setUp(self):
        self.client.force_login(self.user)

    def _get_grid(self, node):
        kwargs = {
            "analysis_id": self.analysis.pk,
            "analysis_version": self.analysis.version,
            "node_id": node.pk,
            "node_version": node.version,
            "extra_filters": "default",
        }
        return self.client.get(reverse("node_data_grid", kwargs=kwargs))

    def test_large_node_shows_banner(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=50_000)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.context["grid_sorting_disabled"])
        self.assertContains(response, "Sorting disabled (50,000 variants)")

    def test_small_node_no_banner(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=500)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.context["grid_sorting_disabled"])
        self.assertNotContains(response, "Sorting disabled")
