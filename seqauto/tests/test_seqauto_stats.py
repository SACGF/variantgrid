"""Enrichment kit stacked bar grouping - how many kits are named and how many fall into "other"."""
import pandas as pd
from django.test import TestCase

from seqauto.seqauto_stats import (
    VARIANTS_TYPE_COL,
    VARIANTS_TYPE_SERIES_COL,
    add_variants_type_series,
    group_enrichment_kits_df,
)


def _kits_df(num_kits: int) -> pd.DataFrame:
    """ One row per sample, kit_00 the biggest so the sort order is predictable """
    rows = []
    for i in range(num_kits):
        rows.extend([{"enrichment_kit__name": f"kit_{i:02d}", "year": 24, "year_month": 2401,
                      "month_offset": 0}] * (num_kits - i))
    return pd.DataFrame(rows)


class EnrichmentKitGroupsTest(TestCase):
    def test_collapses_smallest_kits_into_other(self):
        groups = group_enrichment_kits_df(_kits_df(15), "year", max_groups=10)
        names = [name for name, _ in groups.data]
        self.assertEqual(["kit_00", "kit_01", "kit_02", "kit_03", "kit_04",
                          "kit_05", "kit_06", "kit_07", "kit_08", "other"], names)
        self.assertEqual(6, groups.collapsed_count)
        self.assertEqual([21], dict(groups.data)["other"])  # 6 + 5 + 4 + 3 + 2 + 1
        self.assertIn("top 9", groups.collapsed_help)
        self.assertIn("remaining 6", groups.collapsed_help)

    def test_under_the_limit_has_no_other(self):
        groups = group_enrichment_kits_df(_kits_df(10), "year", max_groups=10)
        self.assertEqual(10, len(groups.data))
        self.assertEqual(0, groups.collapsed_count)
        self.assertIsNone(groups.collapsed_help)


class VariantsTypeSeriesTest(TestCase):
    def test_somatic_types_combine_and_missing_kit_is_unknown(self):
        df = pd.DataFrame({VARIANTS_TYPE_COL: ["G", "M", "S", "U", None],
                           "year": [24] * 5, "year_month": [2401] * 5, "month_offset": [0] * 5})
        add_variants_type_series(df)
        self.assertEqual(["Germline", "Somatic", "Somatic", "Unknown", "Unknown"],
                         list(df[VARIANTS_TYPE_SERIES_COL]))

        groups = group_enrichment_kits_df(df, "year", group_column=VARIANTS_TYPE_SERIES_COL)
        self.assertEqual({"Germline": [1], "Somatic": [2], "Unknown": [2]}, dict(groups.data))

    def test_fake_csv_without_the_column_is_all_unknown(self):
        df = add_variants_type_series(pd.DataFrame({"enrichment_kit__name": ["a", "b"]}))
        self.assertEqual(["Unknown", "Unknown"], list(df[VARIANTS_TYPE_SERIES_COL]))
