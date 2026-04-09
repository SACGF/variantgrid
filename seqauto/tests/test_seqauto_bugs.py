"""
Regression tests for confirmed bugs in the seqauto app.
"""
import os
import shutil
import tempfile
import textwrap

from django.test import SimpleTestCase

from seqauto.illumina.samplesheet import convert_sheet_to_df
from seqauto.models import SequencingRun


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


class TestCheckBasecallsDir(SimpleTestCase):
    """Regression: check_basecalls_dir() always returned False due to missing os.path.join."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_detects_lane_directory(self):
        """Should return True when an L00* subdir exists inside BaseCalls/."""
        basecalls = os.path.join(self.tmpdir, "Data", "Intensities", "BaseCalls")
        os.makedirs(os.path.join(basecalls, "L001"))
        run = SequencingRun.__new__(SequencingRun)
        run.path = self.tmpdir
        self.assertTrue(run.check_basecalls_dir())
