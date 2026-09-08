"""
Copy number on the analysis grid (#1558 §4) - the VCF field binding, and reading the value back out
of the stored CohortGenotype JSON at query time rather than from a packed column.
"""
from django.test import SimpleTestCase

from analysis.grids import VariantGrid
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from snpdb.grid_columns.grid_sample_columns import (
    COPY_NUMBER_COLUMN,
    get_available_format_columns,
    get_copy_number_alias,
    get_variantgrid_zygosity_annotation_kwargs,
)
from snpdb.models import VCFFormat
from snpdb.models.models_enums import VCFInfoTypes
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.views.datatable_view import datatable_response
from upload.vcf.vcf_import import get_copy_number_field


class CopyNumberFieldBindingTest(SimpleTestCase):
    """ Which key the caller wrote copy number under - @see VCF.copy_number_field """

    def test_format_fields_bind_in_preference_order(self):
        self.assertEqual("CN", get_copy_number_field({"CN", "SM", "FC", "GT"}, set(), single_sample=True))
        self.assertEqual("SM", get_copy_number_field({"SM", "FC", "GT"}, set(), single_sample=True))
        self.assertEqual("FC", get_copy_number_field({"FC", "GT"}, set(), single_sample=True))

    def test_info_is_the_fallback_for_a_single_sample_vcf(self):
        """ The Pisces TSO 500 shape - INFO/CN with no FORMAT copy number """
        self.assertEqual("CN", get_copy_number_field({"GT"}, {"CN", "END"}, single_sample=True))

    def test_info_is_not_read_for_a_multi_sample_vcf(self):
        """ One INFO value says nothing about which sample it belongs to """
        self.assertIsNone(get_copy_number_field({"GT"}, {"CN"}, single_sample=False))

    def test_no_copy_number_in_the_header(self):
        self.assertIsNone(get_copy_number_field({"GT", "AD"}, {"END"}, single_sample=True))


class CopyNumberColumnTest(GridExportTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.vcf.copy_number_field = "CN"
        cls.vcf.save()
        VCFFormat.objects.create(vcf=cls.vcf, identifier="CN", number="1",
                                 data_type=VCFInfoTypes.FLOAT, description="Copy number")

        # The proband is index 0 in the VCF's own cohort, which is the FORMAT list's order too
        cgc = cls.cohort.cohort_genotype_collection
        cls.cn_integer_variant = slowly_create_test_variant("1", 3000, "A", "T", cls.genome_build)
        cls._add_genotype(cgc, cls.cn_integer_variant, sample_format=[{"CN": [3]}, {}, {}])
        cls.cn_ratio_variant = slowly_create_test_variant("1", 4000, "A", "T", cls.genome_build)
        cls._add_genotype(cgc, cls.cn_ratio_variant, sample_format=[{"CN": [1.4949999]}, {}, {}])
        cls.no_cn_variant = cls.variants[0]

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(FakeRequest(user=self.user), self.node)
        self.column_name = f"sample_{self.sample.pk}_{COPY_NUMBER_COLUMN}"

    def _columns_by_name(self) -> dict:
        return {rc.name: rc for rc in self.grid.enabled_columns}

    def test_column_is_offered_when_the_vcf_declares_a_field(self):
        self.assertTrue(get_available_format_columns(self.grid.cohorts)[COPY_NUMBER_COLUMN])
        columns = self._columns_by_name()
        self.assertIn(self.column_name, columns)
        # Drawn inside the zygosity cell, so it rides along hidden
        self.assertFalse(columns[self.column_name].visible)
        self.assertEqual("CN proband", columns[self.column_name].label,
                         "labelled with the VCF's own key, not a generic one")

    def test_zygosity_cell_is_told_what_the_key_means(self):
        """ Neither the key nor the header's description of it is on the row """
        zygosity = self._columns_by_name()[f"sample_{self.sample.pk}_samples_zygosity"]
        self.assertEqual({"label": "CN", "title": "Copy number"},
                         zygosity.client_renderer_kwargs["copyNumber"])
        self.assertIn({"label": "Copy number", "column": self.column_name}, zygosity.sort_menu)

    def test_one_annotation_alias_per_sample(self):
        """ There is no packed array to index into, so each sample's value is annotated on its own """
        cgc = self.cohort.cohort_genotype_collection
        annotation_kwargs = get_variantgrid_zygosity_annotation_kwargs([self.cohort], common_variants=True)
        aliases = [get_copy_number_alias(cgc, sample.pk) for sample in self.cohort.get_samples()]
        self.assertTrue(aliases)
        for alias in aliases:
            self.assertIn(alias, annotation_kwargs)
        self.assertEqual([get_copy_number_alias(cgc, self.sample.pk)],
                         self._columns_by_name()[self.column_name].extra_columns)

    def _rows_by_variant_id(self) -> dict:
        self.grid.request.GET = {"length": "100"}
        return {row["id"]: row for row in datatable_response(self.grid)["data"]}

    def test_rows_carry_the_copy_number(self):
        rows = self._rows_by_variant_id()
        self.assertEqual(3, rows[self.cn_integer_variant.pk][self.column_name])
        self.assertIsNone(rows[self.no_cn_variant.pk][self.column_name],
                          "a row the caller wrote no copy number for")

    def test_a_ratio_is_rounded_server_side(self):
        """ Rounding in the cell would leave the CSV saying something else """
        rows = self._rows_by_variant_id()
        self.assertEqual(1.495, rows[self.cn_ratio_variant.pk][self.column_name])

    def test_sorting_is_numeric(self):
        index = next(i for i, rc in enumerate(self.grid.enabled_columns) if rc.name == self.column_name)
        self.grid.request.GET = {"order[0][column]": str(index), "order[0][dir]": "desc"}
        qs = self.grid.ordering(self.grid.get_initial_queryset())
        self.assertEqual([self.cn_integer_variant.pk, self.cn_ratio_variant.pk],
                         [v.pk for v in qs][:2], "3 sorts above 1.495, not below it as text")

    def test_csv_export_includes_it(self):
        header, rows = self._export_csv(self.node)
        index = header.index("CN proband")
        by_variant = {int(row[0]): row for row in rows}
        self.assertEqual("3", by_variant[self.cn_integer_variant.pk][index])
        self.assertEqual("1.495", by_variant[self.cn_ratio_variant.pk][index])
        self.assertEqual("", by_variant[self.no_cn_variant.pk][index])

    def test_vcf_export_round_trips_it_under_the_vcfs_own_key(self):
        lines = self._export_lines(self.node, export_type="vcf")
        self.assertIn('##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">', lines)
        records = [line.split("\t") for line in lines if not line.startswith("#")]
        self.assertEqual("CN", records[0][8].split(":")[-1], "appended after the standard keys")
        calls = {(r[0], int(r[1])): r[9] for r in records}

        def copy_number(variant):
            locus = variant.locus
            return calls[(locus.contig.name, locus.position)].split(":")[-1]

        self.assertEqual("3", copy_number(self.cn_integer_variant))
        self.assertEqual("1.495", copy_number(self.cn_ratio_variant))
        self.assertEqual(".", copy_number(self.no_cn_variant))


class CopyNumberNotDeclaredTest(GridExportTestCase):
    """ The VCF has no copy number field, so nothing about the grid changes """

    def test_column_is_absent(self):
        grid = VariantGrid(FakeRequest(user=self.user), self._sample_node())
        self.assertFalse(get_available_format_columns(grid.cohorts)[COPY_NUMBER_COLUMN])
        names = {rc.name for rc in grid.enabled_columns}
        self.assertNotIn(f"sample_{self.sample.pk}_{COPY_NUMBER_COLUMN}", names)

    def test_vcf_export_declares_no_extra_format(self):
        node = self._sample_node()
        lines = self._export_lines(node, export_type="vcf")
        records = [line.split("\t") for line in lines if not line.startswith("#")]
        self.assertEqual("GT:AD:AF:PL:DP:GQ", records[0][8])
