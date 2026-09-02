from django.test import TestCase

from snpdb.grid_columns.composite_examples import (
    COMPOSITE_EXAMPLE_ROWS,
    SAMPLE_EXAMPLE_PREFIX,
    SAMPLE_EXAMPLE_ROW,
)
from snpdb.models import CohortGenotype, VariantGridColumn


class CompositeExamplesTest(TestCase):
    """ The annotation descriptions page draws a real cell from every example row, so an example that
        has drifted from the catalogue draws a cell missing the member it stopped carrying """

    @classmethod
    def setUpTestData(cls):
        cls.composites = [vgc for vgc in VariantGridColumn.objects.prefetch_related("composite_members__column")
                          if vgc.is_composite]

    def test_example_row_matches_its_members(self):
        for composite in self.composites:
            with self.subTest(composite=composite.pk):
                example = COMPOSITE_EXAMPLE_ROWS.get(composite.pk)
                self.assertIsNotNone(example, f"No example row for composite '{composite.pk}'")
                # The composite's own path is allowed - the Variant cell's value is the variant id
                expected = {m.column.variant_column for m in composite.composite_members.all()}
                allowed = expected | {composite.variant_column}
                self.assertFalse(set(example) - allowed, "example row has keys that aren't members")
                self.assertFalse(expected - set(example), "example row is missing members")

    def test_headline_values_are_set(self):
        """ A blank headline draws an empty cell, which documents nothing """
        for composite in self.composites:
            with self.subTest(composite=composite.pk):
                headline = composite.headline_column.variant_column
                value = COMPOSITE_EXAMPLE_ROWS[composite.pk].get(headline)
                self.assertNotIn(value, (None, ""), f"'{composite.pk}' headline {headline} is blank")

    def test_sample_example_row_keys(self):
        for key in SAMPLE_EXAMPLE_ROW:
            self.assertTrue(key.startswith(SAMPLE_EXAMPLE_PREFIX), key)
            format_column = key[len(SAMPLE_EXAMPLE_PREFIX):]
            self.assertIn(format_column, CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE)
