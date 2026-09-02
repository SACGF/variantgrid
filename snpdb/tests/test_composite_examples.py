from django.db.models import fields
from django.test import TestCase

from library.django_utils import resolve_field_path
from snpdb.grid_columns.composite_examples import (
    COMPOSITE_EXAMPLE_ROWS,
    SAMPLE_EXAMPLE_PREFIX,
    SAMPLE_EXAMPLE_ROW,
)
from snpdb.grids import variant_grid_client_extra
from snpdb.models import CohortGenotype, GenomeBuild, Variant, VariantGridColumn


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

    def test_choice_members_hold_their_display_label(self):
        """ A choice field goes to the client through the choices renderer, so a chip that keys off
            one (the ClinVar stars, the AMP tier) only matches an example holding the label """
        for composite in self.composites:
            example = COMPOSITE_EXAMPLE_ROWS[composite.pk]
            for member in composite.member_columns:
                if not member.model_field:
                    continue
                field = resolve_field_path(Variant._meta, member.variant_column)
                if not (isinstance(field, fields.CharField) and field.choices):
                    continue
                value = example[member.variant_column]
                with self.subTest(column=member.variant_column):
                    self.assertIn(value, [None] + [label for _code, label in field.choices])

    def test_clinvar_stars_lookup_finds_the_example_review_statuses(self):
        """ The stars on the ClinVar chips come out of ctx.extra, keyed the way the column arrives """
        stars = variant_grid_client_extra(GenomeBuild.grch38())["clinvarStars"]
        example = COMPOSITE_EXAMPLE_ROWS["classifications"]
        self.assertEqual(stars[example["clinvar__review_status"]], 2)
        self.assertEqual(stars[example["clinvar__somatic_review_status"]], 1)

    def test_sample_example_row_keys(self):
        for key in SAMPLE_EXAMPLE_ROW:
            self.assertTrue(key.startswith(SAMPLE_EXAMPLE_PREFIX), key)
            format_column = key[len(SAMPLE_EXAMPLE_PREFIX):]
            self.assertIn(format_column, CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE)
