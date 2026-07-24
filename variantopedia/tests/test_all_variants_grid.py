import json

from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase
from django.urls import reverse

from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants
from annotation.models import AnnotationRangeLock, AnnotationRun, VariantAnnotationVersion
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneSymbol, GeneSymbolAlias, GeneSymbolAliasSource
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import AllVariantsFilter, GenomeBuild, Sequence, Variant, VariantZygosityCountCollection
from snpdb.variant_filters import VariantType, get_all_variant_types, get_default_all_variants_filters
from variantopedia.grids import AllVariantsGrid


class AllVariantsGridFilterTest(TestCase):
    """ The All Variants page composes its query from the filter buttons, and refuses to scan the whole
        variant table when nothing selective is chosen - @see issue #1663 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_all_variants_grid")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)

        cls.variant = Variant.objects.get(locus__contig__name='13', locus__position=95839002,
                                          locus__ref__seq='C', alt__seq='T')
        cls.contig = cls.variant.locus.contig
        cls.other_variant = Variant.objects.get(locus__contig__name='15', locus__position=32928050,
                                                locus__ref__seq='C', alt__seq='T')

        # A gene overlapping cls.variant, reachable via an alias of its symbol
        transcript_version = create_fake_transcript_version(cls.genome_build)
        cls.gene = transcript_version.gene_version.gene
        cls.gene_symbol = transcript_version.gene_version.gene_symbol
        cls._create_variant_gene_overlap(cls.variant, cls.gene)

        cls.alias_symbol = GeneSymbol.objects.get_or_create(symbol="FAKEALIAS")[0]
        GeneSymbolAlias.objects.get_or_create(alias=cls.gene_symbol.symbol, gene_symbol=cls.alias_symbol,
                                              source=GeneSymbolAliasSource.NCBI)

        reference_alt = Sequence.objects.get_or_create(seq=Variant.REFERENCE_ALT)[0]
        cls.reference_variant = Variant.objects.create(locus=cls.variant.locus, alt=reference_alt,
                                                       end=cls.variant.end)

    @classmethod
    def _create_variant_gene_overlap(cls, variant, gene):
        """ Annotation partitions are Postgres table inheritance with no insert routing, so an ORM create()
            lands in the parent table, where the partition scoped queries can't see it """
        variant_annotation_version = cls.annotation_version.variant_annotation_version
        annotation_range_lock = AnnotationRangeLock.objects.create(version=variant_annotation_version,
                                                                   min_variant=variant, max_variant=variant, count=1)
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock)
        table_name = variant_annotation_version.get_partition_table(
            base_table_name=VariantAnnotationVersion.VARIANT_GENE_OVERLAP)
        with connection.cursor() as cursor:
            cursor.execute(f'INSERT INTO "{table_name}" (version_id, annotation_run_id, variant_id, gene_id) '
                           'VALUES (%s, %s, %s, %s)',
                           [variant_annotation_version.pk, annotation_run.pk, variant.pk, gene.pk])

    def _grid_variant_ids(self, extra_filters) -> set[int]:
        grid = AllVariantsGrid(self.user, self.genome_build.name, extra_filters=extra_filters)
        qs = grid._get_base_queryset()
        # AbstractVariantGrid.get_queryset() always installs this join, so min_count has its aliases
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        if q := grid._get_q():
            qs = qs.filter(q)
        return set(qs.values_list("pk", flat=True))

    def test_contig_filter(self):
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk]})
        self.assertIn(self.variant.pk, variant_ids)
        self.assertNotIn(self.other_variant.pk, variant_ids)

    def test_unselective_filters_return_nothing(self):
        """ Variant type / min count alone don't restrict the scan, so the grid matches nothing """
        variant_ids = self._grid_variant_ids({"variant_types": [VariantType.SNV], "min_count": 1})
        self.assertEqual(set(), variant_ids)

    def test_no_extra_filters_uses_saved_filters(self):
        """ A direct grid hit (CSV export, bookmarked URL) behaves like the page """
        AllVariantsFilter.objects.update_or_create(user=self.user, genome_build=self.genome_build,
                                                   defaults={"filters": {"contig_ids": [self.contig.pk]}})
        variant_ids = self._grid_variant_ids({})
        self.assertIn(self.variant.pk, variant_ids)
        self.assertNotIn(self.other_variant.pk, variant_ids)

    def test_reference_variants_excluded(self):
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk]})
        self.assertIn(self.variant.pk, variant_ids)
        self.assertNotIn(self.reference_variant.pk, variant_ids)

    def test_variant_type_filter(self):
        """ Selecting only symbolic types drops the SNVs """
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk], "variant_types": ["<DEL>"]})
        self.assertNotIn(self.variant.pk, variant_ids)

    def test_symbolic_type_requires_svlen(self):
        """ A '<DEL>' alt only counts as a symbolic deletion when it has an SVLEN """
        alt = Sequence.objects.get_or_create(seq=VCFSymbolicAllele.DEL)[0]
        locus = self.variant.locus
        symbolic_variant = Variant.objects.create(locus=locus, alt=alt, end=locus.position + 500, svlen=500)
        no_svlen_variant = Variant.objects.create(locus=locus, alt=alt, end=locus.position, svlen=None)

        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk], "variant_types": ["<DEL>"]})
        self.assertIn(symbolic_variant.pk, variant_ids)
        self.assertNotIn(no_svlen_variant.pk, variant_ids)

    def test_all_variant_types_selected_is_no_restriction(self):
        all_types = self._grid_variant_ids({"contig_ids": [self.contig.pk],
                                            "variant_types": get_all_variant_types()})
        no_types = self._grid_variant_ids({"contig_ids": [self.contig.pk]})
        self.assertEqual(no_types, all_types)

    def test_gene_symbol_filter(self):
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk],
                                              "gene_symbols": [self.gene_symbol.symbol]})
        self.assertEqual({self.variant.pk}, variant_ids)

    def test_gene_symbol_filter_traverses_aliases(self):
        """ Searching the alias finds the variants of the symbol it aliases """
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk],
                                              "gene_symbols": [self.alias_symbol.symbol]})
        self.assertEqual({self.variant.pk}, variant_ids)

    def test_gene_symbol_alone_is_selective(self):
        variant_ids = self._grid_variant_ids({"gene_symbols": [self.gene_symbol.symbol]})
        self.assertEqual({self.variant.pk}, variant_ids)

    def test_min_count_filter(self):
        """ min_count >= 1 restricts to variants observed in samples - the test data has none """
        variant_ids = self._grid_variant_ids({"contig_ids": [self.contig.pk], "min_count": 1})
        self.assertEqual(set(), variant_ids)

    def test_default_filters_pick_chr21(self):
        filters = get_default_all_variants_filters(self.genome_build)
        contig_21 = self.genome_build.standard_contigs.get(name="21")
        self.assertEqual([contig_21.pk], filters["contig_ids"])


class AllVariantsGridSortTest(TestCase):
    """ Only the allowlisted columns are sortable - sorting on a joined/unindexed column full-sorts the whole
        result set before LIMIT, blowing the statement_timeout (issues #1279, #1651) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_all_variants_grid_sort")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)

    def _grid(self) -> AllVariantsGrid:
        return AllVariantsGrid(self.user, self.genome_build.name)

    def test_non_allowlisted_sort_falls_back_to_pk(self):
        grid = self._grid()
        sorted_qs = grid._sort_items(grid._get_base_queryset(), sidx="variantannotation__gene_symbol", sord="asc")
        self.assertEqual(["-pk"], list(sorted_qs.query.order_by))

    def test_allowlisted_sort_kept_with_pk_tiebreaker(self):
        grid = self._grid()
        sorted_qs = grid._sort_items(grid._get_base_queryset(), sidx="id", sord="asc")
        order_by = list(sorted_qs.query.order_by)
        self.assertEqual("-pk", order_by[-1])
        self.assertGreater(len(order_by), 1)

    def test_colmodels_outside_allowlist_not_sortable(self):
        grid = self._grid()
        colmodels = grid.get_colmodels()
        self.assertTrue(colmodels)
        for cm in colmodels:
            if cm["name"] in AllVariantsGrid.SORTABLE_FIELDS:
                self.assertNotEqual(False, cm.get("sortable"))
            else:
                self.assertIs(False, cm.get("sortable"))

    def test_default_sort_is_id(self):
        grid = self._grid()
        self.assertEqual("id", grid.extra_config["sortname"])


class AllVariantsFilterPersistenceTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_all_variants_filter_persist")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)

    def setUp(self):
        self.client.force_login(self.user)

    def _post(self, filters):
        url = reverse("set_all_variants_filter", kwargs={"genome_build_name": self.genome_build.name})
        return self.client.post(url, data=json.dumps(filters), content_type="application/json")

    def test_round_trip(self):
        filters = {
            "contig_ids": [1, 2],
            "non_standard_contigs": True,
            "gene_symbols": ["BRCA1"],
            "variant_types": [VariantType.SNV],
            "min_count": 3,
        }
        response = self._post(filters)
        self.assertEqual(200, response.status_code)
        all_variants_filter = AllVariantsFilter.get(self.user, self.genome_build)
        self.assertEqual(filters, all_variants_filter.filters)

    def test_overwrites_previous(self):
        self._post({"contig_ids": [1]})
        self._post({"contig_ids": [2]})
        all_variants_filter = AllVariantsFilter.get(self.user, self.genome_build)
        self.assertEqual([2], all_variants_filter.filters["contig_ids"])

    def test_empty_filters_falls_back_to_defaults(self):
        self._post({"contig_ids": [1]})
        self._post({})
        all_variants_filter = AllVariantsFilter.get(self.user, self.genome_build)
        self.assertEqual({}, all_variants_filter.filters)

    def test_get_requires_post(self):
        url = reverse("set_all_variants_filter", kwargs={"genome_build_name": self.genome_build.name})
        self.assertEqual(405, self.client.get(url).status_code)
