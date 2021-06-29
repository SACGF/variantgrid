from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, CohortNode, CohortNodeZygosityFiltersCollection, FilterNode, FilterNodeItem, \
    SampleNode, GeneListNode, GeneListNodeGeneList
from analysis.tasks.analysis_update_tasks import create_and_launch_analysis_tasks
from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.gene_matching import GeneMatcher
from genes.models import GeneList, GeneListGeneSymbol
from snpdb.models import GenomeBuild, ImportStatus
from snpdb.tests.test_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestCloneAnalysisNodes(TestCase):
    """ It's only worth testing AnalysisNodes with related objects, as normal fields are just copied OK """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        cls.transcript_version = create_fake_transcript_version(cls.grch37)
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.trio = create_fake_trio(cls.user, cls.grch37)

        user = User.objects.get_or_create(username='testuser')[0]
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        sample = cls.trio.get_samples()[0]
        cls.sample_node = SampleNode.objects.create(analysis=cls.analysis, sample=sample)

        # Need to do this so sample node is ready, so children attached will work OK
        create_and_launch_analysis_tasks(cls.analysis.pk, run_async=False)

        # Gene List
        cls.gene_list = GeneList.objects.get_or_create(name="fake list",
                                                   user=cls.analysis.user,
                                                   import_status=ImportStatus.SUCCESS)[0]
        GeneListGeneSymbol.objects.get_or_create(gene_list=cls.gene_list, gene_symbol=cls.gene_symbol)

        gene_annotation_release = cls.annotation_version_grch37.gene_annotation_version.gene_annotation_release
        gm = GeneMatcher(gene_annotation_release)
        gm.match_unmatched_in_hgnc_and_gene_lists()


    def _test_clone_node(self, node):
        clone = node.save_clone()
        self.assertEqual(node.get_q(), clone.get_q(), f"Clone of {node} has same get_q()")

    def _test_make_child_of_sample_node_and_clone(self, node):
        # Make child of sample node
        node.add_parent(self.sample_node)
        node.save()
        # Make clone also child of sample node
        clone = node.save_clone()
        clone.add_parent(self.sample_node)
        clone.save()

        msg = f"SampleNode child node of type '{node.get_class_name()}' clone has same get_q()"
        # Test as strings
        node_q = str(node._get_node_q())
        clone_q = str(clone._get_node_q())
        print(f"Type '{node.get_class_name()}' Node: {node_q} <=> {clone_q}")
        self.assertEqual(node_q, clone_q, msg)

    def test_clone_cohort_node(self):
        cohort = self.trio.cohort
        cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=cohort,
                                                accordion_panel=CohortNode.PER_SAMPLE_ZYGOSITY)

        fc = CohortNodeZygosityFiltersCollection.get_for_node_and_cohort(cohort_node, cohort)
        zf1 = fc.cohortnodezygosityfilter_set.order_by("pk").first()
        zf1.zygosity_hom = False
        zf1.save()

        zf2 = fc.cohortnodezygosityfilter_set.order_by("pk").last()
        zf2.zygosity_het = False
        zf2.save()

        self._test_clone_node(cohort_node)

    def test_clone_filter_node(self):
        filter_node = FilterNode.objects.create(analysis=self.analysis)
        FilterNodeItem.objects.create(filter_node=filter_node, sort_order=1,
                                      operation="eq", field="id", data="42")

        self._test_make_child_of_sample_node_and_clone(filter_node)

    def test_clone_gene_list_node(self):
        gene_list_node = GeneListNode.objects.create(analysis=self.analysis)
        gene_list_node.genelistnodegenelist_set.create(gene_list_node=gene_list_node, gene_list=self.gene_list)
        # TODO: This isn't actually working properly
        self._test_make_child_of_sample_node_and_clone(gene_list_node)


    def test_clone_intersection_node(self):
        pass

    def test_clone_pedigree_node(self):
        pass

    def test_clone_phenotype_node(self):
        pass

    def test_clone_population_node(self):
        pass

    def test_clone_sample_node(self):
        pass

    def test_clone_tag_node(self):
        pass

    def test_clone_trio_node(self):
        pass

    def test_clone_zygosity_node(self):
        pass
