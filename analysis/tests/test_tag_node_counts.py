"""
Tests for per-tag node counts (issue #21) - every tag used in an analysis can have its own count
badge on the nodes, added automatically as variants are tagged - and the TagNode editor's tag
count toggles (issue #1820), which filter the grid on any of the selected tags.
"""
from django.urls import reverse

from analysis.grids import VariantGrid
from analysis.models import Analysis, NodeStatus, TagNode, VariantTag
from analysis.models.enums import TagNodeInput, TagNodeMode, ZygosityNodeZygosity
from analysis.models.nodes.analysis_node import NodeCount
from analysis.models.nodes.filters.zygosity_node import ZygosityNode
from analysis.models.nodes.node_counts import (
    get_extra_filters_count,
    get_extra_filters_q,
    get_node_counts_mine_and_available,
    get_node_extra_filters_q,
)
from analysis.models.nodes.node_utils import update_analysis_tag_node_counts
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from snpdb.models import Tag, Variant
from snpdb.models.models_enums import BuiltInFilters, TagFilter
from snpdb.models.models_user_settings import AbstractNodeCountSettings
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


class TagNodeCountTestCase(GridExportTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.tag = Tag.objects.get_or_create(pk="artefact")[0]
        cls.tag_label = TagFilter.label(cls.tag.pk)
        cls.other_tag = Tag.objects.get_or_create(pk="review")[0]
        cls.other_tag_label = TagFilter.label(cls.other_tag.pk)

    def _tag_variant(self, variant, analysis=None, tag=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=self.genome_build, analysis=analysis or self.analysis,
                                         variant=variant, tag=tag or self.tag, user=self.user)

    def _node_count_types(self) -> list[str]:
        return [label for label, _ in self.analysis.get_node_count_types()]

    def _tag_node(self, mode=TagNodeMode.THIS_ANALYSIS) -> TagNode:
        """ A TagNode over the sample node, holding every tagged variant in its scope """
        parent = self._sample_node()
        node = TagNode.objects.create(analysis=self.analysis, mode=mode)
        node.add_parent(parent)
        node._cached_parents = None  # Clear stale cache from create()'s save()
        node.save()
        node.update(count=node.get_queryset().count(), status=NodeStatus.READY)
        return node


class TestTagNodeCountConfig(TagNodeCountTestCase):
    def test_tagging_adds_the_tag_node_count(self):
        self._tag_variant(self.variants[0])
        self.assertIn(self.tag_label, self._node_count_types())

    def test_untagging_the_last_variant_removes_the_tag_node_count(self):
        variant_tag = self._tag_variant(self.variants[0])
        variant_tag.delete()
        self.assertNotIn(self.tag_label, self._node_count_types())

    def test_untagging_keeps_the_node_count_while_other_variants_are_tagged(self):
        variant_tag = self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1])
        variant_tag.delete()
        self.assertIn(self.tag_label, self._node_count_types())

    def test_auto_add_off_leaves_node_counts_alone(self):
        self.analysis.node_count_auto_add_tags = False
        self.analysis.save()
        self._tag_variant(self.variants[0])
        self.assertNotIn(self.tag_label, self._node_count_types())

    def test_tag_is_offered_in_the_available_tags_column(self):
        _my, available_filters, available_tags = get_node_counts_mine_and_available(self.analysis)
        self.assertNotIn(self.tag_label, {nc["pk"] for nc in available_filters})
        available_by_pk = {nc["pk"]: nc for nc in available_tags}
        self.assertIn(self.tag_label, available_by_pk)
        self.assertEqual(self.tag.pk, available_by_pk[self.tag_label]["description"])
        self.assertEqual(self.tag.pk, available_by_pk[self.tag_label]["tag_id"])

    def test_deleted_tag_drops_out_of_the_node_count_types(self):
        self.analysis.set_node_count_types([BuiltInFilters.TOTAL, TagFilter.label("never_existed")])
        self.assertEqual([BuiltInFilters.TOTAL], self._node_count_types())


class TestTagNodeCountValues(TagNodeCountTestCase):
    """ The recount runs in a task after tagging, so call it directly rather than through celery """

    def _tag_count(self, node) -> int:
        update_analysis_tag_node_counts(self.analysis)
        return NodeCount.objects.get(node_version=node.node_version, label=self.tag_label).count

    def test_counts_tagged_variants_in_the_node(self):
        node = self._sample_node()
        self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1])
        self.assertEqual(2, self._tag_count(node))

    def test_untagging_lowers_the_count(self):
        node = self._sample_node()
        variant_tag = self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1])
        variant_tag.delete()
        self.assertEqual(1, self._tag_count(node))

    def test_variant_outside_the_node_is_not_counted(self):
        node = self._sample_node()
        # No genotype for this one, so the sample node doesn't contain it
        outside_variant = slowly_create_test_variant("3", 12345, "A", "T", self.genome_build)
        self._tag_variant(outside_variant)
        self.assertEqual(0, self._tag_count(node))

    def test_node_whose_parent_is_not_ready_is_skipped(self):
        """ Tagging can land while the analysis is reloading - the child counts these when it loads """
        parent = self._sample_node()
        child = ZygosityNode.objects.create(analysis=self.analysis, sample=self.sample,
                                            zygosity=ZygosityNodeZygosity.HET)
        child.add_parent(parent)
        child._cached_parents = None  # Clear stale cache from create()'s save()
        child.save()
        child.update(status=NodeStatus.READY)
        parent.update(status=NodeStatus.DIRTY)

        self._tag_variant(self.variants[0])
        update_analysis_tag_node_counts(self.analysis)

        self.assertFalse(NodeCount.objects.filter(node_version__node=child, label=self.tag_label).exists())

    def test_recount_bumps_modified_so_the_client_sees_it(self):
        """ nodes_status hands the client counts_modified - it's how it knows a recount landed """
        node = self._sample_node()
        self._tag_variant(self.variants[0])
        update_analysis_tag_node_counts(self.analysis)
        node_count = NodeCount.objects.get(node_version=node.node_version, label=self.tag_label)
        first_modified = node_count.modified

        update_analysis_tag_node_counts(self.analysis)
        node_count.refresh_from_db()
        self.assertGreater(node_count.modified, first_modified)


class TestTagExtraFiltersQ(TagNodeCountTestCase):
    def _filtered_variant_pks(self) -> set:
        q = get_extra_filters_q(self.analysis, self.tag_label)
        return set(Variant.objects.filter(q).values_list("pk", flat=True))

    def test_filters_to_variants_tagged_in_this_analysis(self):
        self._tag_variant(self.variants[0])
        self.assertEqual({self.variants[0].pk}, self._filtered_variant_pks())

    def test_ignores_tags_made_in_another_analysis(self):
        other_analysis = Analysis(genome_build=self.genome_build)
        other_analysis.set_defaults_and_save(self.user)
        self._tag_variant(self.variants[0], analysis=other_analysis)
        self.assertEqual(set(), self._filtered_variant_pks())


class TestNodeCountsSettingsTab(TagNodeCountTestCase):
    def test_saving_stores_the_tag_count_and_the_auto_add_setting(self):
        self.client.force_login(self.user)
        url = reverse("analysis_settings_node_counts_tab", kwargs={"analysis_id": self.analysis.pk})
        # The checkbox is absent from the POST when unticked
        response = self.client.post(url, {"node_counts": f"{BuiltInFilters.TOTAL},{self.tag_label}"})

        self.assertEqual(200, response.status_code)
        self.analysis.refresh_from_db()  # The view saved against its own instance
        self.assertEqual([BuiltInFilters.TOTAL, self.tag_label], self._node_count_types())
        self.assertFalse(self.analysis.node_count_auto_add_tags)


class TestNodeCountLabels(TagNodeCountTestCase):
    def test_tag_node_count_type_carries_its_tag(self):
        types = dict(AbstractNodeCountSettings.get_types_from_labels([self.tag_label]))
        self.assertEqual(self.tag.pk, types[self.tag_label]["label"])
        self.assertEqual(self.tag.pk, types[self.tag_label]["tag"])
        self.assertTrue(types[self.tag_label]["link"])


class TestTagFilterLabels(TagNodeCountTestCase):
    def test_multi_tag_label_round_trips(self):
        label = TagFilter.label_for_tags([self.tag.pk, self.other_tag.pk])
        self.assertEqual([self.tag.pk, self.other_tag.pk], TagFilter.get_tag_ids(label))

    def test_built_in_filter_is_not_a_tag_label(self):
        self.assertEqual([], TagFilter.get_tag_ids(BuiltInFilters.CLINVAR))

    def test_a_non_tag_among_the_parts_makes_it_not_a_tag_label(self):
        self.assertEqual([], TagFilter.get_tag_ids(f"{self.tag_label},{BuiltInFilters.CLINVAR}"))

    def test_multi_tag_label_describes_itself_as_both_tags(self):
        label = TagFilter.label_for_tags([self.tag.pk, self.other_tag.pk])
        description = AbstractNodeCountSettings.get_node_count_description(label)
        self.assertEqual(f"{self.tag.pk}, {self.other_tag.pk}", description)


class TestMultiTagExtraFilters(TagNodeCountTestCase):
    """ A tag counts summary selection goes down as a comma joined label @see TagFilter """

    def setUp(self):
        super().setUp()
        self.both_tags_label = TagFilter.label_for_tags([self.tag.pk, self.other_tag.pk])

    def test_multi_tag_filter_is_the_union_of_its_tags(self):
        self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1], tag=self.other_tag)
        q = get_extra_filters_q(self.analysis, self.both_tags_label)
        self.assertEqual({self.variants[0].pk, self.variants[1].pk},
                         set(Variant.objects.filter(q).values_list("pk", flat=True)))

    def test_grid_filters_and_counts_without_a_stored_node_count(self):
        """ A compound selection never has a NodeCount row - the grid counts it exactly instead """
        node = self._sample_node()
        self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1], tag=self.other_tag)

        grid = VariantGrid(FakeRequest(user=self.user), node, self.both_tags_label)
        self.assertEqual(2, grid._grid_row_count())
        self.assertEqual({self.variants[0].pk, self.variants[1].pk},
                         set(node.get_queryset().filter(grid._get_q()).values_list("pk", flat=True)))

    def test_create_extra_filter_child_makes_a_tag_node_carrying_both_tags(self):
        node = self._sample_node()
        self.client.force_login(self.user)
        url = reverse("create_extra_filter_child", kwargs={"analysis_id": self.analysis.pk,
                                                           "node_id": node.pk,
                                                           "extra_filters": self.both_tags_label})
        response = self.client.post(url)

        self.assertEqual(200, response.status_code)
        tag_node = TagNode.objects.get(analysis=self.analysis)
        self.assertEqual([self.tag.pk, self.other_tag.pk], sorted(tag_node.tag_ids))

    def test_create_extra_filter_child_keeps_a_global_parents_scope(self):
        parent = self._tag_node(mode=TagNodeMode.ALL_TAGS)
        parent.tagged_within_days = 7
        parent.save()
        self.client.force_login(self.user)
        url = reverse("create_extra_filter_child", kwargs={"analysis_id": self.analysis.pk,
                                                           "node_id": parent.pk,
                                                           "extra_filters": self.tag_label})
        response = self.client.post(url)

        self.assertEqual(200, response.status_code)
        child = TagNode.objects.exclude(pk=parent.pk).get(analysis=self.analysis)
        self.assertEqual(TagNodeMode.ALL_TAGS, child.mode)
        self.assertEqual(7, child.tagged_within_days)


class TestGlobalTagNodeCounts(TagNodeCountTestCase):
    """ A global TagNode's tags reach outside the analysis, so it counts and filters them itself """

    def setUp(self):
        super().setUp()
        self.other_analysis = Analysis(genome_build=self.genome_build)
        self.other_analysis.set_defaults_and_save(self.user)

    def _tag_in_other_analysis(self, variant) -> VariantTag:
        """ Tags from another analysis are matched through the allele, as they may be another build """
        allele = create_mock_allele(variant, self.genome_build)
        variant_tag = self._tag_variant(variant, analysis=self.other_analysis)
        variant_tag.allele = allele
        variant_tag.save()
        return variant_tag

    def test_counts_a_tag_made_in_another_analysis(self):
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node(mode=TagNodeMode.ALL_TAGS)
        self.assertEqual([(self.tag.pk, 1)], node.get_tag_counts())

    def test_ignores_the_tagged_within_days_cutoff(self):
        """ The cutoff decides which variants enter the node, not which tags to count """
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node(mode=TagNodeMode.ALL_TAGS)
        node.tagged_within_days = 0
        node.save()
        self.assertEqual([(self.tag.pk, 1)], node.get_tag_counts())

    def test_local_mode_does_not_count_another_analysis_tag(self):
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node()
        self.assertEqual([], node.get_tag_counts())

    def test_local_node_filters_to_this_analysis_only(self):
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node()
        q = get_node_extra_filters_q(node, self.tag_label)
        self.assertEqual(set(), set(node.get_queryset().filter(q).values_list("pk", flat=True)))

    def test_global_node_filters_to_the_tags_it_counted(self):
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node(mode=TagNodeMode.ALL_TAGS)
        q = get_node_extra_filters_q(node, self.tag_label)
        self.assertEqual({self.variants[0].pk},
                         set(node.get_queryset().filter(q).values_list("pk", flat=True)))

    def test_global_node_ignores_its_analysis_scoped_stored_count(self):
        """ NodeCounts are analysis-scoped, so they disagree with a global node's own filter """
        self._tag_in_other_analysis(self.variants[0])
        node = self._tag_node(mode=TagNodeMode.ALL_TAGS)
        NodeCount.objects.create(node_version=node.node_version, label=self.tag_label, count=99)
        self.assertEqual(1, get_extra_filters_count(node, self.tag_label))

    def test_no_extra_filter_shows_the_whole_node(self):
        node = self._tag_node()
        self.assertIsNone(get_extra_filters_count(node, "default"))
        self.assertIsNone(get_node_extra_filters_q(node, "default"))
        grid = VariantGrid(FakeRequest(user=self.user), node, "default")
        self.assertEqual(node.count, grid._grid_row_count())


class TestTagNodeEditorCounts(TagNodeCountTestCase):
    """ The pills above the TagNode form - toggling one reloads the editor at that extra_filters """

    def _editor_html(self, node, extra_filters="default") -> str:
        self.client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                           "analysis_version": self.analysis.version,
                                           "node_id": node.pk,
                                           "node_version": node.version,
                                           "extra_filters": extra_filters})
        response = self.client.get(url)
        self.assertEqual(200, response.status_code)
        return response.content.decode()

    def test_local_mode_renders_a_pill_per_tag_in_the_analysis(self):
        self._tag_variant(self.variants[0])
        node = self._tag_node()
        html = self._editor_html(node)
        self.assertIn(f'data-tag="{self.tag.pk}"', html)
        self.assertNotIn(f'data-tag="{self.other_tag.pk}"', html)

    def test_counts_come_from_the_node_not_the_configured_node_counts(self):
        """ Tags applied before the analysis had that tag's node count still get a pill (#1820) """
        self._tag_variant(self.variants[0])
        self.analysis.set_node_count_types([BuiltInFilters.TOTAL])
        node = self._tag_node()
        html = self._editor_html(node)
        self.assertIn(f'data-tag="{self.tag.pk}"', html)
        self.assertIn('<span class="count">1</span>', html)

    def test_clear_is_disabled_rather_than_hidden_with_nothing_selected(self):
        """ A hidden button would shift the pills along as tags are toggled """
        self._tag_variant(self.variants[0])
        node = self._tag_node()
        self.assertIn('clear-tags" disabled', self._editor_html(node))
        self.assertNotIn('clear-tags" disabled', self._editor_html(node, self.tag_label))

    def test_selected_tags_come_back_selected(self):
        self._tag_variant(self.variants[0])
        self._tag_variant(self.variants[1], tag=self.other_tag)
        node = self._tag_node()
        html = self._editor_html(node, TagFilter.label_for_tags([self.tag.pk, self.other_tag.pk]))
        self.assertEqual(2, html.count("summary-count tagged-"))
        self.assertEqual(2, html.count(" selected\""))

    def test_excluding_tagged_variants_has_no_counts(self):
        self._tag_variant(self.variants[0])
        node = self._tag_node()
        node.node_input = TagNodeInput.PARENT_NOT_TAGGED
        node.save()
        self.assertNotIn("tag-counts-summary", self._editor_html(node))
