from django.test import SimpleTestCase

from library.utils.text_utils import delimited_row, spreadsheet_safe


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

    def test_delimited_row_quotes_formula_text(self):
        """ Classification and other CSV downloads are written through delimited_row (SACGF/variantgrid_private#3946) """
        row = delimited_row(["=1+1", "+SUM(A1)", "-1.5", "c.-12A>G", 7, None], include_new_line=False)
        self.assertEqual("'=1+1,'+SUM(A1),-1.5,c.-12A>G,7,", row)
