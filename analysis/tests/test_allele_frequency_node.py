from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, NodeAlleleFrequencyFilter
from analysis.models.nodes.filters.allele_frequency_node import AlleleFrequencyNode
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_data import get_fake_annotation_version
from snpdb.fake_data import create_fake_trio, make_cohort_genotype
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


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
        trio = create_fake_trio(user, cls.grch37)
        cls.sample = trio.get_samples()[0]
        cls.cgc = trio.cohort.cohort_genotype_collection

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

    def test_no_sample_or_patient_is_a_configuration_error(self):
        node = AlleleFrequencyNode.objects.create(analysis=self.analysis)
        self.assertIn("No sample or patient selected.", node._get_configuration_errors())

    def test_restricted_range_filters_on_the_samples_allele_frequency(self):
        """ Only the proband's AF is read - the parents' values would pass or fail the other way """
        variants = {}
        for i, (name, proband_af) in enumerate([("in_range", 0.5), ("below", 0.1), ("above", 0.9)]):
            variant = slowly_create_test_variant("1", 1000 + i * 100, "A", "T", self.grch37)
            parents_af = 0.1 if proband_af == 0.5 else 0.5
            make_cohort_genotype(self.cgc, variant, "EEE", allele_frequency=[proband_af, parents_af, parents_af])
            variants[variant.pk] = name

        parent = SampleNode.objects.create(analysis=self.analysis, sample=self.sample)
        status, count = parent.node_counts()
        parent.update(status=status, count=count)
        node = self._node()
        node.add_parent(parent)
        node._cached_parents = None  # Clear stale cache from create()'s save()
        node.save()
        self._set_range(node, 0.15, 0.8)

        filtered_pks = node.get_queryset().filter(pk__in=variants).values_list("pk", flat=True)
        self.assertEqual({"in_range"}, {variants[pk] for pk in filtered_pks})
