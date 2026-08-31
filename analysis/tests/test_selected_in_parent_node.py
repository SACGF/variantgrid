from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode
from analysis.models.nodes.filters.selected_in_parent_node import NodeVariant, SelectedInParentNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import create_fake_variants
from snpdb.models import Variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestSelectedInParentNode(AnalysisSetupMixin, TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        create_fake_variants(cls.grch37)
        cls.variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])

    def test_filters_to_variants_selected_in_parent(self):
        parent = AllVariantsNode.objects.create(analysis=self.analysis)
        node = SelectedInParentNode.objects.create(analysis=self.analysis)
        node.add_parent(parent)
        node.save()
        node = SelectedInParentNode.objects.get(pk=node.pk)  # pick up the new edge

        selected = self.variants[:2]
        for variant in selected:
            NodeVariant.objects.create(node=parent.analysisnode_ptr, variant=variant)

        matched = set(Variant.objects.filter(node._get_node_q())
                      .filter(pk__in=[v.pk for v in self.variants])
                      .values_list("pk", flat=True))
        self.assertEqual(matched, {v.pk for v in selected})

    def test_always_filters(self):
        """ Nothing selected means nothing passes - the parent is never returned unmodified """
        node = SelectedInParentNode(analysis=self.analysis)
        self.assertTrue(node.modifies_parents())
