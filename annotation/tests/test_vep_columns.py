from django.test import TestCase

from annotation.models.models_enums import (
    ColumnAnnotationCategory,
    VariantAnnotationPipelineType,
    VEPCustom,
    VEPPlugin,
)
from annotation.vep_columns import (
    VEP_COLUMNS,
    _aloft,
    _dbnsfp_v2,
    _gnomad4,
    all_variant_grid_column_ids,
    filter_for,
)
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

    def test_gnomad4_helper_defaults(self):
        c = _gnomad4('AF_afr', 'gnomad_afr_af')
        self.assertEqual(c.source_field, 'AF_afr')
        self.assertEqual(c.variant_grid_columns, ('gnomad_afr_af',))
        self.assertEqual(c.category, ColumnAnnotationCategory.FREQUENCY_DATA)
        self.assertEqual(c.vep_custom, VEPCustom.GNOMAD_4)
        self.assertTrue(c.source_field_has_custom_prefix)
        self.assertEqual(c.genome_builds, frozenset({'GRCh38', 'T2T-CHM13v2.0'}))
        self.assertEqual(c.pipeline_types, frozenset({VariantAnnotationPipelineType.STANDARD}))
        self.assertEqual(c.min_columns_version, 3)

    def test_gnomad4_helper_override(self):
        """ The faf95/faf99 rows narrow to GRCh38 only; helper override must win. """
        c = _gnomad4('faf95', 'gnomad_faf95', genome_builds=frozenset({'GRCh38'}))
        self.assertEqual(c.genome_builds, frozenset({'GRCh38'}))
        self.assertEqual(c.min_columns_version, 3)
        self.assertEqual(c.vep_custom, VEPCustom.GNOMAD_4)

    def test_aloft_carries_canonical_description(self):
        c = _aloft('Aloft_pred', 'aloft_pred')
        self.assertEqual(c.vep_plugin, VEPPlugin.DBNSFP)
        self.assertEqual(c.min_columns_version, 2)
        self.assertTrue(
            c.source_field_processing_description.startswith(
                "Most damaging transcript prediction chosen"
            )
        )

    def test_dbnsfp_v2_helper_defaults(self):
        c = _dbnsfp_v2('REVEL_rankscore', 'revel_rankscore')
        self.assertEqual(c.vep_plugin, VEPPlugin.DBNSFP)
        self.assertEqual(c.category, ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS)
        self.assertEqual(c.genome_builds, frozenset({'GRCh37', 'GRCh38'}))
        self.assertEqual(c.min_columns_version, 2)

    def test_filter_for_gnomad4_grch38_v3_count(self):
        rows = filter_for(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_4,
            gnomad4_minor_version="4.1",
        )
        self.assertEqual(len(rows), 26)

    def test_filter_for_gnomad3_excluded_at_v3(self):
        rows = filter_for(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_3,
        )
        self.assertEqual(rows, ())
