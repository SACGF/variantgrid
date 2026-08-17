from django.test import TestCase

from library.genomics.vcf_utils import UnsortedVCFError, VCFSortChecker


def _check_all(records):
    """ records: (contig, position) - chrom/line_number only matter for the error message """
    sort_checker = VCFSortChecker()
    for line_number, (contig, position) in enumerate(records, start=1):
        sort_checker.check(contig, position, str(contig), line_number)


class TestVCFSortChecker(TestCase):
    def test_sorted(self):
        _check_all([(1, 100), (1, 100), (1, 200), (2, 50), (2, 60), (3, 10)])

    def test_position_goes_backwards(self):
        with self.assertRaises(UnsortedVCFError):
            _check_all([(1, 154560359), (1, 154560345)])

    def test_contig_not_grouped(self):
        with self.assertRaises(UnsortedVCFError):
            _check_all([(1, 100), (2, 100), (1, 200)])

    def test_contig_order_does_not_matter(self):
        """ bcftools index only needs contigs contiguous, not in header order """
        _check_all([(3, 100), (1, 100), (2, 100)])

    def test_error_message_has_positions(self):
        with self.assertRaises(UnsortedVCFError) as cm:
            _check_all([(1, 200), (1, 100)])
        message = str(cm.exception)
        self.assertIn("1:100", message)
        self.assertIn("1:200", message)
        self.assertIn("Line 2", message)
