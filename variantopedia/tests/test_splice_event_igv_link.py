"""The junction a gene-level splice variant offers an IGV link at, having no coordinate of its own (#1908)."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.fake_annotation import get_fake_annotation_version
from genes.models import HGNC, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_splice_event_variant
from snpdb.models import GenomeBuild

AR_HGNC_ID = 644


class TestSpliceEventIgvLink(TestCase):
    """ The page hands the junction's breakpoints to the client, which decides whether IGV links
        show at all - @see renderIgvLocusLinks in grid.js """

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
        self.assertContains(response, 'data-locus="X:66905968-66914514"')

    def test_grid_row_detail_offers_the_junction(self):
        url = reverse("variant_grid_row_detail",
                      kwargs={"variant_id": self.variant.pk,
                              "annotation_version_id": self.annotation_version.pk})
        self.assertContains(self.client.get(url), 'data-locus="X:66905968-66914514"')
