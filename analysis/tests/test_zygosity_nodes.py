from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, CohortNode, AnalysisNode, AllVariantsNode
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestZygosityNodes(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        cls.trio = create_fake_trio(user, cls.grch37)

        samples_list = list(cls.trio.cohort.get_samples()[:2])  # Leave out last sample to be sub-cohort
        cls.sub_cohort = cls.trio.cohort.create_sub_cohort(user, samples_list)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def test_cohort_node_zyg_filters(self):
        cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=self.trio.cohort,
                                                accordion_panel=CohortNode.COUNT)
        self._test_zygosity(cohort_node)

    def test_sub_cohort_node_zyg_filters(self):
        sub_cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=self.sub_cohort,
                                                accordion_panel=CohortNode.COUNT)
        self._test_zygosity(sub_cohort_node)

    def test_all_variants_node_zyg_filters(self):
        all_variants_node = AllVariantsNode(analysis=self.analysis)
        self._test_zygosity(all_variants_node)

    def _test_zygosity(self, node: AnalysisNode):
        node.min_het_or_hom_count = 1
        node.max_het_or_hom_count = 99
        node.min_ref_count = 1
        node.max_ref_count = 99
        node.min_het_count = 1
        node.max_het_count = 99
        node.min_hom_count = 1
        node.max_hom_count = 99
        node.save()  # Also makes a whole bunch of needed stuff

        qs = node.get_queryset()
        count = qs.count()
        self.assertIsInstance(count, int)
