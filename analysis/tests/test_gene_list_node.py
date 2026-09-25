from datetime import timedelta

from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.models import GeneListNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import create_fake_variants
from annotation.models import AnnotationRun, VariantGeneOverlap
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol, PanelAppPanel, PanelAppPanelLocalCache, PanelAppServer
from pathtests.models import PathologyTest, PathologyTestVersion
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

    def _panel_app_node(self, cached_version: str) -> GeneListNode:
        server = PanelAppServer.objects.create(name="Test PanelApp", url="https://panelapp.example.com",
                                               icon_css_class="")
        panel = PanelAppPanel.objects.create(server=server, panel_id=1, disease_group="", disease_sub_group="",
                                             name="Test Panel", status="public", current_version="0.2")
        local_cache = PanelAppPanelLocalCache.objects.create(panel_app_panel=panel, version=cached_version)
        node = GeneListNode.objects.create(analysis=self.analysis, accordion_panel=GeneListNode.PANEL_APP_GENE_LIST)
        node.genelistnodepanelapppanel_set.create(panel_app_panel=panel, panel_app_panel_local_cache=local_cache)
        return node

    def test_panel_app_older_cached_version_warns(self):
        warnings = self._panel_app_node(cached_version="0.1").get_warnings()
        self.assertTrue(any("v.0.1 while latest is 0.2" in w for w in warnings))

    def test_panel_app_current_cached_version_no_warning(self):
        self.assertEqual([], self._panel_app_node(cached_version="0.2").get_warnings())

    @override_settings(PATHOLOGY_TEST_STALE_WARNING_DAYS=365)
    def test_stale_pathology_test_warns(self):
        pathology_test = PathologyTest.objects.create(name="test")
        ptv = PathologyTestVersion.objects.create(pathology_test=pathology_test, gene_list=self.gene_list)
        node = GeneListNode.objects.create(analysis=self.analysis, pathology_test_version=ptv,
                                           accordion_panel=GeneListNode.PATHOLOGY_TEST_GENE_LIST)
        self.assertFalse(any("last modified" in w for w in node.get_warnings()))

        PathologyTestVersion.objects.filter(pk=ptv.pk).update(modified=timezone.now() - timedelta(days=400))
        node = GeneListNode.objects.get(pk=node.pk)
        self.assertTrue(any("last modified 400 days ago" in w for w in node.get_warnings()))
