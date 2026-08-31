from django.test import TestCase, override_settings

from analysis.models import GeneListNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import create_fake_variants
from annotation.models import AnnotationRun, VariantGeneOverlap
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol
from snpdb.models import ImportStatus, Variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestGeneListNode(AnalysisSetupMixin, TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        annotation_version = cls.analysis.annotation_version
        release = annotation_version.gene_annotation_version.gene_annotation_release
        transcript_version = create_fake_transcript_version(cls.grch37, release=release)

        cls.gene_list = GeneList.objects.create(name="fake list", user=cls.analysis.user,
                                                import_status=ImportStatus.SUCCESS)
        GeneListGeneSymbol.objects.create(gene_list=cls.gene_list,
                                          gene_symbol=transcript_version.gene_version.gene_symbol)

        create_fake_variants(cls.grch37)
        cls.variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        cls.in_gene_list, cls.other_1, cls.other_2 = cls.variants
        VariantGeneOverlap.objects.create(version=annotation_version.variant_annotation_version,
                                          annotation_run=AnnotationRun.objects.create(),
                                          gene=transcript_version.gene_version.gene,
                                          variant=cls.in_gene_list)

    def _node(self, **kwargs) -> GeneListNode:
        node = GeneListNode.objects.create(analysis=self.analysis, **kwargs)
        node.genelistnodegenelist_set.create(gene_list=self.gene_list)
        return node

    def _matched(self, node) -> set:
        return set(Variant.objects.filter(node._get_node_q())
                   .filter(pk__in=[v.pk for v in self.variants])
                   .values_list("pk", flat=True))

    def test_no_gene_list_does_not_modify_parents(self):
        self.assertFalse(GeneListNode.objects.create(analysis=self.analysis).modifies_parents())

    def test_gene_list_modifies_parents(self):
        self.assertTrue(self._node().modifies_parents())

    def test_filters_to_variants_overlapping_gene_list(self):
        self.assertEqual(self._matched(self._node()), {self.in_gene_list.pk})

    def test_exclude_inverts_the_filter(self):
        self.assertEqual(self._matched(self._node(exclude=True)), {self.other_1.pk, self.other_2.pk})

    def test_exclude_has_no_known_contigs(self):
        """ Excluding a gene list can match anything, so contig optimisation has to be skipped """
        self.assertIsNone(self._node(exclude=True)._get_node_contigs())
