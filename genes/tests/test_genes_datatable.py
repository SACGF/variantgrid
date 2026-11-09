from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse

from annotation.fake_data import get_fake_annotation_version
from annotation.models import GeneAnnotation
from genes.enums import AnnotationConsortium
from genes.fake_data import create_fake_gene_version, create_fake_transcript_version
from genes.grids import GenesColumns
from genes.models import ReleaseGeneVersion, TranscriptVersion
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import GenomeBuild


class GenesDatatableTests(TestCase):
    """ The genes page filters by release and by "column is/isn't null" shortcuts """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="genes_datatable_user")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.transcript_version = create_fake_transcript_version(cls.genome_build)

    def _config(self, **params) -> GenesColumns:
        url = reverse('genes_datatable', kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, params)
        request.resolver_match = resolve(url)
        request.user = self.user
        return GenesColumns(request)

    def test_is_null_column_filter(self):
        annotation_version = get_fake_annotation_version(self.genome_build)
        gene_annotation_version = annotation_version.gene_annotation_version
        for gene_id, symbol, omim_terms in [("ENSG90000000001", "WITH_OMIM", "OMIM:123456"),
                                            ("ENSG90000000002", "NO_OMIM", None)]:
            gene_version = create_fake_gene_version(self.genome_build, gene_id, symbol, AnnotationConsortium.ENSEMBL)
            ReleaseGeneVersion.objects.create(release=gene_annotation_version.gene_annotation_release,
                                              gene_version=gene_version)
            with temporary_db_table(GeneAnnotation, gene_annotation_version.get_partition_table()):
                GeneAnnotation.objects.create(version=gene_annotation_version, gene=gene_version.gene,
                                              omim_terms=omim_terms)

        for is_null, expected_symbol in [("true", "NO_OMIM"), ("false", "WITH_OMIM")]:
            with self.subTest(is_null=is_null):
                config = self._config(column="gene_version__gene__geneannotation__omim_terms", is_null=is_null)
                qs = config.filter_queryset(config.get_initial_queryset())
                self.assertEqual({expected_symbol}, set(qs.values_list("gene_version__gene_symbol", flat=True)))

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

    def test_gene_version_ids_in_region_boundaries(self):
        """ cdot exons are 0-based half-open, regions (cytobands) 1-based inclusive. The fake transcript's
            exons span [34787800, 35049344), ie 1-based 34787801-35049344 """
        tv = self.transcript_version
        for start, end, overlaps in [(1, 34787800, False), (1, 34787801, True),
                                     (35049344, 40000000, True), (35049345, 40000000, False),
                                     (34900000, 34900001, True)]:
            with self.subTest(start=start, end=end):
                gene_version_ids = TranscriptVersion.get_gene_version_ids_in_region(self.genome_build, tv.contig,
                                                                                    start, end)
                self.assertEqual(overlaps, tv.gene_version_id in set(gene_version_ids))

    def test_contig_region_filter(self):
        """ region_start/region_end narrow a contig's genes to a band ('start' is the DataTables offset) """
        tv = self.transcript_version
        release = get_fake_annotation_version(self.genome_build).gene_annotation_version.gene_annotation_release
        ReleaseGeneVersion.objects.get_or_create(release=release, gene_version=tv.gene_version)
        for region, expected in [({}, {"RUNX1"}),
                                 ({"region_start": 34800000, "region_end": 34800100}, {"RUNX1"}),
                                 ({"region_start": 1, "region_end": 1000}, set())]:
            with self.subTest(region=region):
                config = self._config(contig_id=tv.contig_id, **region)
                qs = config.filter_queryset(ReleaseGeneVersion.objects.filter(release=release))
                self.assertEqual(expected, set(qs.values_list("gene_version__gene_symbol", flat=True)))
