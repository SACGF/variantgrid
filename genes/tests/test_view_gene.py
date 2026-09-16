"""view_gene for the genes the template can't link: legacy (fake) genes and versions with no symbol."""
from django.contrib.auth.models import User
from django.urls import reverse

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import Gene, GeneVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils.unittest_utils import URLTestCase, _make_test_client
from snpdb.models.models_genome import GenomeBuild


class ViewGeneTest(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.gene_version = create_fake_transcript_version(cls.grch37).gene_version

    def setUp(self):
        self.client = _make_test_client()
        self.client.force_login(self.user)

    def _get_gene(self, gene_id):
        return self.client.get(reverse("view_gene", kwargs={"gene_id": gene_id}))

    def test_legacy_gene_redirects_to_symbol(self):
        symbol = self.gene_version.gene_symbol_id
        Gene.objects.create(identifier=Gene.FAKE_GENE_ID_PREFIX + symbol,
                            annotation_consortium=AnnotationConsortium.ENSEMBL)
        response = self._get_gene(Gene.FAKE_GENE_ID_PREFIX + symbol)
        self.assertRedirects(response, reverse("view_gene_symbol", kwargs={"gene_symbol": symbol}),
                             fetch_redirect_response=False)

    def test_legacy_gene_without_symbol_is_404(self):
        Gene.objects.create(identifier=Gene.FAKE_GENE_ID_PREFIX, annotation_consortium=AnnotationConsortium.ENSEMBL)
        self.assertEqual(self._get_gene(Gene.FAKE_GENE_ID_PREFIX).status_code, 404)

    def test_version_without_symbol_renders_no_symbol_link(self):
        gene = Gene.objects.create(identifier="ENSG00000238009", annotation_consortium=AnnotationConsortium.ENSEMBL)
        GeneVersion.objects.create(gene=gene, version=1, genome_build=self.grch37,
                                   import_source=self.gene_version.import_source)
        response = self._get_gene(gene.pk)
        self.assertEqual(response.status_code, 200)
        self.assertNotIn(b"view_gene_symbol/None", response.content)
