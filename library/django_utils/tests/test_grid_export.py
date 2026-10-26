from django.test import SimpleTestCase

from library.django_utils.grid_export import spreadsheet_safe


class SpreadsheetSafeTest(SimpleTestCase):
    def test_formula_text_is_quoted(self):
        cases = [
            ('=HYPERLINK("x")', '\'=HYPERLINK("x")'),
            ("@SUM(A1)", "'@SUM(A1)"),
            ("-1.5", "-1.5"),  # a number
            ("-", "-"),  # VEP's empty value
            ("c.-12A>G", "c.-12A>G"),
            (-3, -3),
            (None, None),
        ]
        for value, expected in cases:
            self.assertEqual(expected, spreadsheet_safe(value), value)
