from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, NodeAlleleFrequencyFilter
from analysis.models.nodes.filters.allele_frequency_node import AlleleFrequencyNode
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestAlleleFrequencyNode(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username="test_TestAlleleFrequencyNode")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)
        cls.sample = create_fake_trio(user, cls.grch37).get_samples()[0]

    def _node(self) -> AlleleFrequencyNode:
        return AlleleFrequencyNode.objects.create(analysis=self.analysis, sample=self.sample)

    @staticmethod
    def _set_range(node, min_af, max_af):
        naff = NodeAlleleFrequencyFilter.objects.get(node=node)
        naff.nodeallelefrequencyrange_set.all().delete()
        naff.nodeallelefrequencyrange_set.create(min=min_af, max=max_af)

    def test_new_node_has_an_unrestricted_filter(self):
        node = self._node()
        af_range = NodeAlleleFrequencyFilter.objects.get(node=node).nodeallelefrequencyrange_set.get()
        self.assertEqual((af_range.min, af_range.max), (0, 100))

    def test_unrestricted_range_does_not_filter(self):
        self.assertFalse(self._node().modifies_parents())

    def test_no_sample_is_a_configuration_error(self):
        node = AlleleFrequencyNode.objects.create(analysis=self.analysis)
        self.assertIn("No sample selected.", node._get_configuration_errors())

    def test_restricted_range_filters_on_the_samples_allele_frequency(self):
        node = self._node()
        self._set_range(node, 0.15, 0.8)
        arg_q_dict = node._get_node_arg_q_dict()
        q_strs = [str(q) for q_dict in arg_q_dict.values() for q in q_dict.values()]
        self.assertTrue(any("allele_frequency__0__gte" in s for s in q_strs), q_strs)
        self.assertTrue(any("allele_frequency__0__lte" in s for s in q_strs), q_strs)
