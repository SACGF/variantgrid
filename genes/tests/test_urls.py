import json
import uuid

from django.contrib.auth.models import User
import unittest

from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import CachedWebResource
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import CanonicalTranscriptCollection, GeneCoverageCollection, GeneCoverage, GeneList, PanelAppPanel, \
    GeneListCategory, Gene, CanonicalTranscript, GeneListGeneSymbol, PanelAppServer
from genes.models_enums import AnnotationConsortium
from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from snpdb.models import ImportStatus, DataState
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.test_data import create_fake_trio


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        trio = create_fake_trio(cls.user_owner, cls.grch37)
        cls.cohort = trio.cohort
        cls.transcript_version = create_fake_transcript_version(cls.grch37)
        cls.transcript = cls.transcript_version.transcript
        cls.gene = cls.transcript_version.gene_version.gene
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.uncovered_gene = Gene.objects.get_or_create(identifier="fake_gene",
                                                        annotation_consortium=AnnotationConsortium.ENSEMBL)[0]

        ac = AnnotationConsortium.ENSEMBL
        ctc, _ = CanonicalTranscriptCollection.objects.get_or_create(description="fake canonical transcripts",
                                                                     filename="/tmp/foo.txt",
                                                                     annotation_consortium=ac,
                                                                     file_md5sum="not_a_real_hash",
                                                                     genome_build=cls.grch37)

        canonical_transcript = CanonicalTranscript.objects.get_or_create(collection=ctc,
                                                                         gene_symbol=cls.gene_symbol,
                                                                         transcript=cls.transcript)[0]

        gcg, created = GeneCoverageCollection.objects.get_or_create(path="/tmp/foo4.txt",
                                                                    data_state=DataState.COMPLETE,
                                                                    genome_build=cls.grch37)
        if created:
            GeneCoverage.objects.create(gene_coverage_collection=gcg,
                                        original_gene_symbol="orig_symbol",
                                        original_transcript_id="orig_refseq",
                                        min=41,
                                        mean=50,
                                        std_dev=1.11,
                                        percent_0x=0.0,
                                        percent_10x=100.0,
                                        percent_20x=100.0,
                                        percent_100x=88.2,
                                        sensitivity=41)

            GeneCoverage.objects.create(gene_coverage_collection=gcg,
                                        original_gene_symbol="orig_symbol",
                                        original_transcript_id="orig_refseq",
                                        min=1,
                                        mean=22,
                                        std_dev=1.11,
                                        percent_0x=0.0,
                                        percent_10x=99.0,
                                        percent_20x=89.9,
                                        percent_100x=50.2,
                                        sensitivity=41)

        cls.canonical_transcript_collection = ctc
        cls.gene_coverage_collection = gcg
        cls.gene_list_category = GeneListCategory.objects.get_or_create(name="Fake Category")[0]
        cls.gene_list = GeneList.objects.get_or_create(name="fake list",
                                                       user=cls.user_owner,
                                                       import_status=ImportStatus.SUCCESS)[0]
        GeneListGeneSymbol.objects.create(gene_list=cls.gene_list, gene_symbol=cls.gene_symbol)
        cls.gene_list_w_category = GeneList.objects.get_or_create(name="fake list",
                                                                  category=cls.gene_list_category,
                                                                  user=cls.user_owner,
                                                                  import_status=ImportStatus.SUCCESS)[0]

        _ = CachedWebResource.objects.get_or_create(name="Fake PanelAppPanels", import_status=ImportStatus.SUCCESS)[0]

        server = PanelAppServer.objects.order_by("pk").first()
        cls.panel_app_panel = PanelAppPanel.objects.get_or_create(server=server,
                                                                  panel_id=42,
                                                                  disease_group='Tumour syndromes',
                                                                  disease_sub_group='Tumour syndromes',
                                                                  name=uuid.uuid4(),
                                                                  current_version='1.20')[0]

        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ('view_gene_list', {'gene_list_id': cls.gene_list.pk}, 200),
            ('api_view_gene_list', {"pk": cls.gene_list.pk}, 200),
            ('cohort_hotspot_graph', {"cohort_id": cls.cohort.pk, "transcript_id": cls.transcript.pk}, 200),
        ]

        cls.PRIVATE_AUTOCOMPLETE_URLS = [
            ('category_gene_list_autocomplete', cls.gene_list, {"q": cls.gene_list.name}),
            ('gene_list_autocomplete', cls.gene_list, {"q": cls.gene_list.name}),
        ]

        cls.PRIVATE_GRID_LIST_URLS = [
            ("gene_lists_grid", {"gene_list_category_id": cls.gene_list_category.pk}, cls.gene_list_w_category),
            ("uncategorised_gene_lists_grid", {}, cls.gene_list),
            ("gene_lists_for_gene_symbol_grid", {"gene_symbol": cls.gene_symbol.pk}, cls.gene_list),
            ("gene_list_genes_grid", {"gene_list_id": cls.gene_list.pk}, ("gene_symbol__symbol", cls.gene_symbol)),
            ("canonical_transcript_collections_grid", {}, cls.canonical_transcript_collection),
            ("canonical_transcript_collection_grid", {"pk": cls.canonical_transcript_collection.pk}, ("transcript__identifier", canonical_transcript.transcript)),

            #("gene_coverage_collection_gene_list_grid", coverage_kwargs, cls.gene),
            #("uncovered_genes_grid", coverage_kwargs, cls.uncovered_gene),
        ]

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("genes", {}, 200),
            ("gene_lists_tab", {}, 200),
            ("gene_lists", {}, 200),
            ("qc_coverage", {}, 200),
            ("gene_grid", {}, 200),
            ("canonical_transcripts", {}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testNoPermissionUrls(self):
        gene_symbol_kwargs = {"gene_symbol": self.gene_symbol}
        URL_NAMES_AND_KWARGS = [
            ("genome_build_genes", {"genome_build_name": self.grch37.name}, 200),
            ("view_gene", {"gene_id": self.gene.pk}, 200),
            ("view_gene_symbol", gene_symbol_kwargs, 200),
            ("view_transcript", {"transcript_id": self.transcript.pk}, 200),
            ("view_transcript_version", {"transcript_id": self.transcript_version.transcript_id,
                                         "version": self.transcript_version.version}, 200),
            # Test accession with version
            ("view_transcript_accession", {"transcript_accession": self.transcript_version.accession}, 200),
            # Test accession has no version
            ("view_transcript_accession", {"transcript_accession": self.transcript_version.transcript_id}, 200),
            # ("api_panel_app_gene_evidence", gene_symbol_kwargs, 200),
            ("api_gene_info", gene_symbol_kwargs, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testGridUrls(self):
        """ Grids w/o permissions """
        genome_build_name = self.grch37.name

        GRID_LIST_URLS = [
            ("gene_symbol_variants_grid", {"gene_symbol": self.gene_symbol.pk, "genome_build_name": genome_build_name, "op": "config"}, 200),
            ("genes_grid", {"genome_build_name": genome_build_name, "op": "config"}, 200),
            ("gene_coverage_collection_grid", {"gene_coverage_collection_id": self.gene_coverage_collection.pk,
                                               "op": "config"}, 200),
        ]
        self._test_urls(GRID_LIST_URLS, self.user_non_owner)

    def testAutocompleteUrls(self):
        """ Autocompletes w/o permissions """
        panel_app_forward = json.dumps({"server_id": self.panel_app_panel.server_id})
        AUTOCOMPLETE_URLS = [
            ('panel_app_panel_autocomplete', self.panel_app_panel, {"q": self.panel_app_panel.name,
                                                                    "forward": panel_app_forward}),
            ('gene_autocomplete', self.gene, {"q": self.gene_symbol}),
            ('transcript_autocomplete', self.transcript, {"q": self.transcript.identifier}),
            ('gene_symbol_autocomplete', self.gene_symbol, {"q": self.gene_symbol}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user_non_owner, True)

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)

    def testAutocompletePermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testAutocompleteNoPermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_non_owner, False)

    def testGridListPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testGridListNoPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_non_owner, False)

# TODO: Still to test
# commented out PRIVATE_GRID_LIST_URLS
# ('gene_lists/create_custom', views.create_custom_gene_list, name='create_custom_gene_list'),
# ('view_canonical_transcript_collection/<pk>', views.view_canonical_transcript_collection, name='view_canonical_transcript_collection'),
# ('qc_coverage/<genome_build_name>', views.qc_coverage, name='genome_build_qc_coverage'),
# ('gene_coverage_collection_graphs/<genome_build_name>/<slug:gene_id>', views.gene_coverage_collection_graphs, name='gene_coverage_collection_graphs'),
# ('qc_gene_list_coverage_graphs/<genome_build_name>/<int:gene_list_id>', views.qc_gene_list_coverage_graphs, name='qc_gene_list_coverage_graphs'),
# ('gene_grid/<path:columns_from_url>', views.gene_grid, name='passed_gene_grid'),
# ('hotspot_graph/gene/<genome_build_name>/<gene_id>',
#               views.HotspotGraphView.as_view(), name='gene_hotspot_graph'),
# ('hotspot_graph/transcript/<genome_build_name>/<transcript_id>',
#               views.HotspotGraphView.as_view(), name='transcript_hotspot_graph'),
# ('hotspot_graph/classifications/gene/<genome_build_name>/<gene_id>',
#               views.ClassificationsHotspotGraphView.as_view(), name='classifications_gene_hotspot_graph'),
# ('hotspot_graph/classifications/transcript/<genome_build_name>/<transcript_id>',
#               views.ClassificationsHotspotGraphView.as_view(), name='classifications_transcript_hotspot_graph'),
# ('hotspot_graph/cohort/<int:cohort_id>/<transcript_id>',
#               views.CohortHotspotGraphView.as_view(), name='cohort_hotspot_graph'),
# ('hotspot_graph/public', views.PublicRUNX1HotspotGraphView.as_view(), name='public_hotspot_graph'),
#
# ('api/gene_list/modify/<int:pk>', views_rest.ModifyGeneListView.as_view(), name='api_modify_gene_list'),
# ('api/gene_list/create', views_rest.CreateGeneListView.as_view(), name='api_create_gene_list'),
# ('api/named_gene_list/<category__name>/<name>', views_rest.NamedGeneListView.as_view(), name='api_named_gene_list'),
# ('api/gene/batch_info', views_rest.BatchGeneInfoView.as_view(), name='api_batch_gene_info'),
# ('api/text_to_gene_list', views_rest.TextToGeneListView.as_view(), name='api_text_to_gene_list'),


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
