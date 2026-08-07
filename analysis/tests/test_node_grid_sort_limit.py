from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.grids import VariantGrid
from analysis.models import Analysis, AllVariantsNode
from annotation.fake_annotation import get_fake_annotation_version
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

    def _grid(self, count):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=count)
        return VariantGrid(self.user, node)

    def test_small_node_sorting_enabled(self):
        grid = self._grid(count=500)
        self.assertFalse(grid.sorting_disabled())

    def test_large_node_sorting_disabled(self):
        grid = self._grid(count=10_000)
        self.assertTrue(grid.sorting_disabled())

    def test_unknown_count_sorting_disabled(self):
        grid = self._grid(count=None)
        self.assertTrue(grid.sorting_disabled())

    def test_large_node_ignores_requested_sort(self):
        grid = self._grid(count=50_000)
        qs = grid._get_base_queryset()
        sorted_qs = grid._sort_items(qs, sidx="variantannotation__gene_symbol", sord="asc")
        self.assertEqual(list(sorted_qs.query.order_by), ["-pk"])

    def test_small_node_keeps_requested_sort(self):
        grid = self._grid(count=500)
        qs = grid._get_base_queryset()
        sorted_qs = grid._sort_items(qs, sidx="start", sord="asc")
        # Requested column first, PK tiebreaker last
        self.assertEqual(list(sorted_qs.query.order_by)[-1], "-pk")
        self.assertGreater(len(sorted_qs.query.order_by), 1)

    def test_large_node_colmodels_not_sortable(self):
        grid = self._grid(count=50_000)
        colmodels = grid.get_colmodels()
        self.assertTrue(colmodels)
        self.assertTrue(all(cm.get("sortable") is False for cm in colmodels))

    def test_small_node_colmodels_sortable_default(self):
        grid = self._grid(count=500)
        colmodels = grid.get_colmodels()
        # Sorting isn't force-disabled - columns keep their normal (non-False) sortable setting
        self.assertFalse(all(cm.get("sortable") is False for cm in colmodels))

    def test_large_node_no_initial_sortname(self):
        grid = self._grid(count=50_000)
        self.assertNotIn("sortname", grid.extra_config)


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
        self.assertContains(response, "Sorting is disabled")

    def test_small_node_no_banner(self):
        node = AllVariantsNode.objects.create(analysis=self.analysis, count=500)
        response = self._get_grid(node)
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.context["grid_sorting_disabled"])
        self.assertNotContains(response, "Sorting is disabled")
