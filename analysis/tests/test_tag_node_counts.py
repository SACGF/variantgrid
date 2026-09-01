"""
Tests for per-tag node counts (issue #21) - every tag used in an analysis can have its own count
badge on the nodes, added automatically as variants are tagged.
"""
from django.urls import reverse

from analysis.models import Analysis, NodeStatus, VariantTag
from analysis.models.enums import ZygosityNodeZygosity
from analysis.models.nodes.analysis_node import NodeCount
from analysis.models.nodes.filters.zygosity_node import ZygosityNode
from analysis.models.nodes.node_counts import get_extra_filters_q, get_node_counts_mine_and_available
from analysis.models.nodes.node_utils import update_analysis_tag_node_counts
from analysis.tests.test_grid_export import GridExportTestCase
from snpdb.models import Tag, Variant
from snpdb.models.models_enums import BuiltInFilters, TagFilter
from snpdb.models.models_user_settings import AbstractNodeCountSettings
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class TagNodeCountTestCase(GridExportTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.tag = Tag.objects.get_or_create(pk="artefact")[0]
        cls.tag_label = TagFilter.label(cls.tag.pk)

    def _tag_variant(self, variant, analysis=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=self.genome_build, analysis=analysis or self.analysis,
                                         variant=variant, tag=self.tag, user=self.user)

    def _node_count_types(self) -> list[str]:
        return [label for label, _ in self.analysis.get_node_count_types()]


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

    def test_tag_is_offered_as_an_available_node_count(self):
        _my, available = get_node_counts_mine_and_available(self.analysis)
        available_by_pk = {nc["pk"]: nc for nc in available}
        self.assertIn(self.tag_label, available_by_pk)
        self.assertEqual(self.tag.pk, available_by_pk[self.tag_label]["description"])

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
