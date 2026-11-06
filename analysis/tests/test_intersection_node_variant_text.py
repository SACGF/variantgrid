from unittest.mock import Mock, patch

from django.test import TestCase, override_settings

from analysis.forms.forms_nodes import IntersectionNodeForm
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.tests.utils import AnalysisSetupMixin
from analysis.variant_text import parse_region, resolve_variant_text
from library.preview_request import PreviewData
from snpdb.models.models_variant import Variant
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
            self.search_data = search_data
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

    def test_cytoband_resolved_to_region_without_search(self):
        resolution = self._resolve("11p13\n11p99")
        self.assertEqual(resolution.regions, ["11:31000001-36400000"])
        self.assertEqual(resolution.unresolved, ["11p99"])
        self.search_data.assert_not_called()

    def test_line_matching_nothing_is_reported_not_dropped(self):
        resolution = self._resolve("rs99999999")
        self.assertEqual(resolution.variant_ids, [])
        self.assertEqual(resolution.unresolved, ["rs99999999"])

    def test_dbsnp_ids_searched_one_at_a_time_until_found(self):
        with patch("analysis.variant_text.search_data") as search_data:
            search_data.side_effect = lambda _user, search_string: Mock(
                results=[_search_result(self.variant)] if search_string == "rs2150111024" else [])
            resolution = resolve_variant_text(self.analysis.user, self.grch37, "rs11958868&rs2150111024&rs3")
        self.assertEqual(resolution.variant_ids, [self.variant.pk])
        self.assertEqual([c.args[1] for c in search_data.call_args_list], ["rs11958868", "rs2150111024"])

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

    def test_coordinates_in_the_database_resolved_without_search(self):
        """ IGV, VCF and gnomAD formats are looked up together rather than searched line by line (#2531) """
        resolution = self._resolve("1:169519049 T>C\n1 169519050 . G A\n1-169519049-T-C")
        self.assertEqual(resolution.variant_ids, [self.variant.pk, self.other_variant.pk])
        self.search_data.assert_not_called()

    def test_coordinate_not_in_the_database_goes_to_search(self):
        resolution = self._resolve("1:169519049 T>G")
        self.assertEqual(resolution.unresolved, ["1:169519049 T>G"])
        self.search_data.assert_called_once()

    @override_settings(SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN=True)
    def test_classified_only_search_is_left_to_search(self):
        """ On Shariant a non-admin's search only finds classified variants, so a coordinate in the database
            isn't resolved locally - when search doesn't find it, neither does the node """
        resolution = self._resolve("1:169519049 T>C")
        self.assertEqual(resolution.variant_ids, [])
        self.assertEqual(resolution.unresolved, ["1:169519049 T>C"])


class IntersectionNodeVariantsQTest(AnalysisSetupMixin, TestCase):
    def test_no_entries_matches_nothing(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="rs99999999")
        self.assertEqual(node._get_node_q(), node.q_none())

    def test_variant_ids_and_regions_are_ored(self):
        listed = slowly_create_test_variant("1", 5000, "A", "T", self.grch37)
        in_region = slowly_create_test_variant("1", 10050, "A", "T", self.grch37)
        elsewhere = slowly_create_test_variant("1", 20000, "A", "T", self.grch37)
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="x\ny", variant_ids=[listed.pk], variant_regions=["1:10000-10100"])
        pks = [listed.pk, in_region.pk, elsewhere.pk]
        filtered_pks = set(Variant.objects.filter(pk__in=pks).filter(node._get_node_q()).values_list("pk", flat=True))
        self.assertEqual({listed.pk, in_region.pk}, filtered_pks)

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

    def test_summary_shows_what_a_cytoband_resolved_to(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="11p13", variant_regions=["11:31000001-36400000"])
        self.assertIn("11p13 = 11:31000001-36400000", node._get_variant_text_summary())

    def test_form_accepts_cytoband(self):
        node = IntersectionNode.objects.create(analysis=self.analysis)
        form = IntersectionNodeForm({"variant_text": "11p13", "accordion_panel": IntersectionNode.VARIANTS},
                                    instance=node)
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()
        self.assertEqual(node.variant_regions, ["11:31000001-36400000"])
        self.assertEqual(node.variant_text_unresolved, [])

    def test_warnings_name_the_unresolved_entries(self):
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.VARIANTS,
                                variant_text="rs6025\nrs99999999", variant_ids=[1],
                                variant_text_unresolved=["rs99999999"])
        self.assertTrue(any("rs99999999" in w for w in node.get_warnings()))
