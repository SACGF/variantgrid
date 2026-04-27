from django.test import TestCase

from annotation.models.models_enums import VariantAnnotationPipelineType, VEPCustom
from annotation.vep_columns import VEP_COLUMNS, all_variant_grid_column_ids, filter_for
from snpdb.models import VariantGridColumn


class VepColumnsRegistryTest(TestCase):

    def test_referenced_variant_grid_columns_exist(self):
        known = set(VariantGridColumn.objects.values_list("pk", flat=True))
        referenced = all_variant_grid_column_ids()
        missing = referenced - known
        self.assertFalse(missing, f"vep_columns.py references unknown VariantGridColumn ids: {sorted(missing)}")

    def test_dedup_keys_unique(self):
        seen = set()
        dupes = []
        for c in VEP_COLUMNS:
            k = (c.source_field, c.vep_plugin, c.vep_custom,
                 c.min_columns_version, c.max_columns_version,
                 c.min_vep_version, c.max_vep_version,
                 c.genome_builds, c.pipeline_types,
                 c.variant_grid_columns)
            if k in seen:
                dupes.append(k)
            seen.add(k)
        self.assertFalse(dupes, f"Duplicate VEPColumnDef keys: {dupes}")

    def test_filter_for_returns_subset(self):
        all_count = len(VEP_COLUMNS)
        filtered = filter_for(genome_build_name="GRCh37")
        self.assertLessEqual(len(filtered), all_count)
        self.assertGreater(len(filtered), 0)

    def test_gnomad4_minor_version_filter(self):
        common_kwargs = dict(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_4,
        )
        rows_40 = [c for c in filter_for(gnomad4_minor_version="4.0", **common_kwargs)
                   if "gnomad_filtered" in c.variant_grid_columns]
        rows_41 = [c for c in filter_for(gnomad4_minor_version="4.1", **common_kwargs)
                   if "gnomad_filtered" in c.variant_grid_columns]
        self.assertEqual(len(rows_40), 1)
        self.assertEqual(rows_40[0].source_field, "gnomad_filtered")
        self.assertEqual(len(rows_41), 1)
        self.assertEqual(rows_41[0].source_field, "FILTER")
