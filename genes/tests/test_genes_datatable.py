from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse

from annotation.fake_data import get_fake_annotation_version
from genes.fake_data import create_fake_transcript_version
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

    def test_columns_drawn_as_the_variant_grid_draws_them(self):
        """ A gene field is the gene level of a VariantGridColumn, and takes its label, tooltip and
            the variant grid's cell renderer; the ones with none still get a width for the fixed layout """
        config = self._config()
        columns = {rc.key: rc for rc in config.enabled_columns}
        interactions = columns["gene_version__gene__geneannotation__dbnsfp_gene__interactions_consensus_pathdb"]
        self.assertEqual(interactions.label, "ConsensusPathDB interactions")
        self.assertEqual(interactions.client_renderer, "VariantGridFormat.bracketedTerms")
        self.assertTrue(interactions.header_title)
        pathways = columns["gene_version__gene__geneannotation__dbnsfp_gene__pathway_consensus_pathdb"]
        self.assertEqual(pathways.client_renderer, "VariantGridFormat.separated")
        self.assertEqual(pathways.client_renderer_kwargs, {"separator": ";"})
        self.assertEqual(columns["gene_version__gene_symbol__symbol"].client_renderer, "renderGeneSymbol")
        self.assertEqual(columns["gene_version__version"].label, "Version")
        self.assertTrue(config.table_class)
        self.assertTrue(config.default_column_width)
