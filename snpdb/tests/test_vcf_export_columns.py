from django.test import TestCase

from analysis.grid_export import _get_column_vcf_info
from snpdb.models import VariantGridColumn
from snpdb.models.models_enums import VCFInfoTypes
from snpdb.vcf_export_columns import COLUMN_VCF_INFO


class TestVCFExportColumns(TestCase):
    """ COLUMN_VCF_INFO replaces the old ColumnVCFInfo DB table (see snpdb migration
        0195_delete_columnvcfinfo). These guard the invariants the table used to enforce. """

    def test_info_ids_unique(self):
        info_ids = [c.info_id for c in COLUMN_VCF_INFO]
        self.assertEqual(len(info_ids), len(set(info_ids)))

    def test_columns_unique(self):
        # Preserves the old OneToOne(VariantGridColumn) guarantee - one info row per column
        columns = [c.column for c in COLUMN_VCF_INFO]
        self.assertEqual(len(columns), len(set(columns)))

    def test_types_are_valid(self):
        for c in COLUMN_VCF_INFO:
            self.assertIsInstance(c.type, VCFInfoTypes)

    def test_consumer_keys_by_variant_column_with_label_type(self):
        # The VCF export contract: entries keyed by the column's variant_column (query path),
        # with type rendered as its human label (e.g. 'Float').
        variant_column_by_name = dict(
            VariantGridColumn.objects.values_list("grid_column_name", "variant_column")
        )
        column_vcf_info = _get_column_vcf_info()

        af_1kg_def = next(c for c in COLUMN_VCF_INFO if c.info_id == "1KG_AF")
        af_1kg = column_vcf_info[variant_column_by_name[af_1kg_def.column]]
        self.assertEqual(af_1kg["info_id"], "1KG_AF")
        self.assertEqual(af_1kg["type"], "Float")
        self.assertEqual(af_1kg["number"], 1)
        self.assertEqual(af_1kg["column__variant_column"], variant_column_by_name[af_1kg_def.column])

        gnomad_filtered_def = next(c for c in COLUMN_VCF_INFO if c.info_id == "GNOMAD_FILTERED")
        gnomad_filtered = column_vcf_info[variant_column_by_name[gnomad_filtered_def.column]]
        self.assertEqual(gnomad_filtered["info_id"], "GNOMAD_FILTERED")
        self.assertEqual(gnomad_filtered["type"], "Flag")

    def test_consumer_skips_absent_columns(self):
        # A def whose grid column isn't present on this deployment is simply skipped
        # (matching the old CASCADE-delete behaviour).
        af_1kg_def = next(c for c in COLUMN_VCF_INFO if c.info_id == "1KG_AF")
        variant_column = VariantGridColumn.objects.get(pk=af_1kg_def.column).variant_column
        VariantGridColumn.objects.filter(pk=af_1kg_def.column).delete()

        column_vcf_info = _get_column_vcf_info()
        self.assertNotIn(variant_column, column_vcf_info)
