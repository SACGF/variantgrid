"""
Regression tests for confirmed bugs in the seqauto app.
"""
import os
import shutil
import tempfile
import textwrap

from django.test import SimpleTestCase

from library.utils.date_utils import parse_yymm
from seqauto.illumina.samplesheet import convert_sheet_to_df
from seqauto.seqauto_stats import get_sequencing_run_yymm


def _write_samplesheet(tmpdir, content, run_dir="200920_NB501009_0410_AHNLYFBGXG"):
    run_path = os.path.join(tmpdir, run_dir)
    os.makedirs(run_path, exist_ok=True)
    sheet_path = os.path.join(run_path, "SampleSheet.csv")
    with open(sheet_path, "w") as f:
        f.write(content)
    return sheet_path


class TestSampleSheetNaNBarcode(SimpleTestCase):
    """Regression: partial Index2 column produced 'CAGATC|nan' barcodes."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_partial_index2_does_not_produce_nan_barcodes(self):
        """When Index2 column exists but some rows are empty, barcode must not contain 'nan'."""
        content = textwrap.dedent("""\
            [Header]
            [Data]
            Sample_ID,Sample_Name,Index,Index2
            SAMPLE_A,SampleA,GCCAAT,ATCGAT
            SAMPLE_B,SampleB,CAGATC,
        """)
        path = _write_samplesheet(self.tmpdir, content)
        df = convert_sheet_to_df(path)
        for bc in df["barcode"]:
            self.assertNotIn(
                "nan", str(bc).lower(),
                f"Barcode should not contain 'nan' for missing Index2, got: {bc!r}",
            )


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
