"""
Regression tests for confirmed bugs in the seqauto app.
"""
from django.test import SimpleTestCase

from library.utils.date_utils import parse_yymm
from seqauto.seqauto_stats import get_sequencing_run_yymm


class TestSequencingRunYYMM(SimpleTestCase):
    def test_parses_illumina_run_date(self):
        self.assertEqual(get_sequencing_run_yymm("250730_M01234_0001_000000000-ABCDE"), "2507")

    def test_dash_separator(self):
        self.assertEqual(get_sequencing_run_yymm("250730-M01234_0001_000000000-ABCDE"), "2507")

    def test_skips_non_date_name(self):
        self.assertIsNone(get_sequencing_run_yymm("some_run_without_a_date"))
        self.assertIsNone(get_sequencing_run_yymm(None))


class TestSampleEnrichmentKitsDF(SimpleTestCase):
    def test_skips_runs_with_dodgy_dates(self):
        """ Placeholder names like 000000_M01234_... have leading digits that aren't a real date - the graph
            drops those rows rather than plotting them decades before the rest """
        DODGY = ["000000_M01234_0001_000000000-ABCDE",  # YYMM 0000
                 "000715_M01234_0001_000000000-ABCDE",  # YYMM 0007 - read as July 2000
                 "001234_M01234_0001_000000000-ABCDE",  # YYMM 0012
                 "251330_M01234_0001_000000000-ABCDE"]  # month 13
        for name in DODGY:
            with self.assertRaises(ValueError, msg=name):
                parse_yymm(get_sequencing_run_yymm(name))
