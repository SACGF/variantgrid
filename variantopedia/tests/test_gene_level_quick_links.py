"""The links a gene-level variant offers, having no coordinate of its own: IGV at a splice junction (#1908)
or at the genes of a fusion, and CIViC for a junction CIViC records (#1909)."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.fake_data import get_fake_annotation_version
from genes.models import HGNC, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.tests.gene_level_test_utils import create_splice_event_variant
from snpdb.models import GenomeBuild

AR_HGNC_ID = 644


class TestGeneLevelQuickLinks(TestCase):
    """ The page hands Quick Links what a gene-level variant's own coordinate can't give it, and the
        client decides whether IGV links show at all - @see VCLinks in vc_links.js """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.genome_build = GenomeBuild.grch37()
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        hgnc_import = HGNCImport.objects.create()
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="AR")

    def setUp(self):
        self.client.force_login(self.user)
        self.variant = create_splice_event_variant("AR", "V7").variant

    def test_variant_page_offers_the_junction(self):
        response = self.client.get(reverse("view_variant", kwargs={"variant_id": self.variant.pk}))
        self.assertContains(response, 'data.link_data.igv_locus = "X:66905968-66914514"')

    def test_fusion_page_offers_both_genes(self):
        """ Space separated, which IGV opens split screen """
        fusion = create_gene_fusion("BCR", "ABL1")
        response = self.client.get(reverse("view_variant", kwargs={"variant_id": fusion.variant.pk}))
        self.assertContains(response, 'data.link_data.igv_locus = "BCR ABL1"')

    def test_grid_row_detail_offers_the_junction(self):
        url = reverse("variant_grid_row_detail",
                      kwargs={"variant_id": self.variant.pk,
                              "annotation_version_id": self.annotation_version.pk})
        self.assertContains(self.client.get(url), 'data-locus="X:66905968-66914514"')

    def test_variant_page_links_the_junctions_civic_variant(self):
        response = self.client.get(reverse("view_variant", kwargs={"variant_id": self.variant.pk}))
        self.assertContains(response, 'data.link_data.civic_variant_url = "https://civicdb.org/variants/362"')
