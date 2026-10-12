import os
import tempfile

from django.test import TestCase, override_settings

from snpdb.models import GenomeBuild
from upload.vcf.vcf_ref_check import (
    RefMismatchCount,
    count_ref_mismatches,
    get_ref_mismatch_message,
    read_snvs,
)

VCF_TEXT = """##fileformat=VCFv4.1
##contig=<ID=1,length=249250621>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1\t.\tA\tG\t.\t.\t.
1\t2\t.\tAC\tA\t.\t.\t.
1\t3\t.\tC\tA,T\t.\t.\t.
1\t4\t.\tN\tT\t.\t.\t.
1\t5\t.\tc\tt\t.\t.\t.
"""


@override_settings(VCF_IMPORT_REF_MISMATCH_WARN_FRACTION=0.05, VCF_IMPORT_REF_MISMATCH_FAIL_FRACTION=0.5)
class TestVCFRefCheck(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.grch37()
        cls.grch38 = GenomeBuild.grch38()

    def test_read_snvs_only_biallelic_acgt(self):
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f:
            f.write(VCF_TEXT)
        try:
            self.assertEqual(read_snvs(f.name, 10), [("1", 1, "A"), ("1", 5, "C")])
            self.assertEqual(read_snvs(f.name, 1), [("1", 1, "A")])
        finally:
            os.unlink(f.name)

    def test_count_skips_unknown_contig_and_ambiguous_reference(self):
        fasta = {"1": "ACGTNacgt"}
        snvs = [("1", 1, "A"), ("1", 2, "A"), ("1", 5, "A"), ("1", 6, "A"), ("1", 50, "A"), ("MT", 1, "A")]
        count = count_ref_mismatches(snvs, self.grch37, fasta)
        # pos 5 is N and pos 50 past the end - neither counted; MT not in fasta
        self.assertEqual((count.num_checked, count.num_mismatched), (3, 1))

    def test_message_thresholds(self):
        self.assertIsNone(get_ref_mismatch_message(RefMismatchCount(self.grch37, 100, 5), []))

        message, fail = get_ref_mismatch_message(RefMismatchCount(self.grch37, 100, 10), [])
        self.assertFalse(fail)

        other = [RefMismatchCount(self.grch38, 100, 0)]
        message, fail = get_ref_mismatch_message(RefMismatchCount(self.grch37, 100, 75), other)
        self.assertTrue(fail)
        self.assertIn("75%", message)
        self.assertIn("called against GRCh38", message)
        self.assertIn("genome_build", message)

    def test_too_few_snvs_to_fail(self):
        _, fail = get_ref_mismatch_message(RefMismatchCount(self.grch37, 3, 3), [])
        self.assertFalse(fail)
