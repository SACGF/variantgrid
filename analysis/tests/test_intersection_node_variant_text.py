from unittest.mock import patch

from django.test import TestCase

from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.tests.utils import AnalysisSetupMixin
from analysis.variant_text import parse_region, resolve_variant_text
from library.preview_request import PreviewData
from snpdb.search import SearchResult
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


def _search_result(obj) -> SearchResult:
    return SearchResult(preview=PreviewData.for_object(obj))


class ResolveVariantTextTest(AnalysisSetupMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.variant = slowly_create_test_variant("1", 169519049, "T", "C", cls.grch37)
        cls.other_variant = slowly_create_test_variant("1", 169519050, "G", "A", cls.grch37)

    def _resolve(self, text, search_results=None):
        with patch("analysis.variant_text.search_data") as search_data:
            search_data.return_value.results = search_results or []
            return resolve_variant_text(self.analysis.user, self.grch37, text)

    def test_blank_lines_ignored(self):
        resolution = self._resolve("\n   \n")
        self.assertEqual(resolution.variant_ids, [])
        self.assertEqual(resolution.unresolved, [])

    def test_region_resolved_without_search(self):
        resolution = self._resolve("chr1:169,519,049-169520049")
        self.assertEqual(resolution.regions, ["1:169519049-169520049"])
        self.assertEqual(resolution.unresolved, [])

    def test_region_unknown_contig_is_an_error(self):
        resolution = self._resolve("chrBANANA:100-200")
        self.assertEqual(resolution.regions, [])
        self.assertEqual(len(resolution.errors), 1)

    def test_region_start_after_end_is_an_error(self):
        resolution = self._resolve("chr1:200-100")
        self.assertEqual(resolution.regions, [])
        self.assertEqual(len(resolution.errors), 1)

    def test_line_matching_nothing_is_reported_not_dropped(self):
        resolution = self._resolve("rs99999999")
        self.assertEqual(resolution.variant_ids, [])
        self.assertEqual(resolution.unresolved, ["rs99999999"])

    def test_variants_deduplicated_across_lines(self):
        resolution = self._resolve("rs6025\nrs6025", search_results=[_search_result(self.variant)])
        self.assertEqual(resolution.variant_ids, [self.variant.pk])

    def test_non_variant_search_results_ignored(self):
        """ search_data returns everything the search box would - genes, samples, other builds """
        resolution = self._resolve("BRCA1", search_results=[_search_result(self.analysis)])
        self.assertEqual(resolution.variant_ids, [])
        self.assertEqual(resolution.unresolved, ["BRCA1"])

    def test_parse_region_round_trip(self):
        resolution = self._resolve("1:100-200")
        self.assertEqual(parse_region(resolution.regions[0]), ("1", 100, 200))


class IntersectionNodeVariantsQTest(AnalysisSetupMixin, TestCase):
    def test_no_entries_matches_nothing(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="rs99999999")
        self.assertEqual(node._get_node_q(), node.q_none())

    def test_variant_ids_and_regions_are_ored(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="x\ny", variant_ids=[1, 2], variant_regions=["1:100-200"])
        self.assertEqual(node._get_node_q().connector, "OR")

    def test_text_that_resolved_to_nothing_still_modifies_parents(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="rs99999999", variant_text_unresolved=["rs99999999"])
        self.assertTrue(node.modifies_parents())

    def test_large_list_uses_cache(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="x", variant_ids=list(range(IntersectionNode.VARIANTS_CACHE_MIN + 1)))
        self.assertTrue(node.use_cache)

    def test_small_list_does_not_use_cache(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="x", variant_ids=[1, 2, 3])
        self.assertFalse(node.use_cache)

    def test_warnings_name_the_unresolved_entries(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="rs6025\nrs99999999", variant_ids=[1],
                                variant_text_unresolved=["rs99999999"])
        self.assertTrue(any("rs99999999" in w for w in node.get_warnings()))
