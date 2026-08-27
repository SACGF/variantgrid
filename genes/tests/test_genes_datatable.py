from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.grids import GenesColumns
from snpdb.models import GenomeBuild


class GenesDatatableTests(TestCase):
    """ The genes page filters by release and by "column is/isn't null" shortcuts """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="genes_datatable_user")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_transcript_version(cls.genome_build)

    def _config(self, **params) -> GenesColumns:
        url = reverse('genes_datatable', kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, params)
        request.resolver_match = resolve(url)
        request.user = self.user
        return GenesColumns(request)

    def test_is_null_column_filter(self):
        config = self._config(column="gene_version__gene__geneannotation__omim_terms", is_null="true")
        sql = str(config.filter_queryset(config.get_initial_queryset()).query)
        self.assertIn("omim_terms", sql)

    def test_unknown_column_refused(self):
        config = self._config(column="gene_version__gene__geneannotation__omim_terms; DROP TABLE")
        with self.assertRaises(PermissionDenied):
            config.filter_queryset(config.get_initial_queryset())
