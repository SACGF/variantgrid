from django.test import SimpleTestCase

from library.date_utils import get_month_and_year
from seqauto.seqauto_stats import get_sequencing_run_yymm


class TestSequencingRunYYMM(SimpleTestCase):
    def test_parses_illumina_run_date(self):
        self.assertEqual(get_sequencing_run_yymm("250730_M01234_0001_000000000-ABCDE"), "2507")

    def test_dash_separator(self):
        self.assertEqual(get_sequencing_run_yymm("250730-M01234_0001_000000000-ABCDE"), "2507")

    def test_skips_placeholder_date(self):
        """ Zero/invalid dates aren't dates - they used to parse to a YYMM that blew up get_month_and_year """
        for name in ["000000_M01234_0001_000000000-ABCDE",
                     "001234_M01234_0001_000000000-ABCDE",
                     "251330_M01234_0001_000000000-ABCDE",
                     "250700_M01234_0001_000000000-ABCDE"]:
            self.assertIsNone(get_sequencing_run_yymm(name), name)

    def test_skips_non_date_name(self):
        self.assertIsNone(get_sequencing_run_yymm("some_run_without_a_date"))
        self.assertIsNone(get_sequencing_run_yymm(None))


class TestGetMonthAndYear(SimpleTestCase):
    def test_yymm(self):
        self.assertEqual(get_month_and_year("2201"), (1, 22))

    def test_leading_zero_year(self):
        """ int() drops the leading zero - 0907 must still be July 2009 """
        self.assertEqual(get_month_and_year(907), (7, 9))

    def test_invalid_yyyymm_format_raises(self):
        with self.assertRaises(ValueError):
            get_month_and_year("202201")
