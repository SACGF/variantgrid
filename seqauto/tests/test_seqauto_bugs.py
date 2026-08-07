"""
Regression tests for confirmed bugs in the seqauto app.
"""
import os
import shutil
import tempfile
import textwrap

from django.test import SimpleTestCase

from seqauto.illumina.samplesheet import convert_sheet_to_df


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
