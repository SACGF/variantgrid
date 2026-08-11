"""
Undeclared FILTER values can't just be passed through - bcftools norm dies on one rather than warning
("Invalid BCF, the FILTER tag id=N ... not present in the header"), and vcf_clean_and_filter strips it.
write_cleaned_vcf_header declares them instead, which is what lets the value reach the grid.
"""
from django.test import TestCase

from library.genomics.vcf_utils import scan_undeclared_filter_ids
from snpdb.models import VCFFilter


def _record(filter_column: str) -> str:
    return f"chr1\t100\t.\tT\t<DEL>\t0\t{filter_column}\tSVTYPE=DEL;END=200\tAD:DP\t1:80\n"


class TestScanUndeclaredFilterIds(TestCase):
    def test_finds_undeclared_value_in_compound_filter(self):
        """ SpliceVariants.vcf declares LowQ only, so LowUniqueAlignments is silently deleted today -
            and it's the filter separating the two PASS oncology calls from the background """
        records = [_record("LowQ;LowUniqueAlignments")]
        self.assertEqual(["LowUniqueAlignments"], scan_undeclared_filter_ids(records, {"LowQ"}))

    def test_declared_values_are_left_alone(self):
        records = [_record("LowQ"), _record("PASS"), _record(".")]
        self.assertEqual([], scan_undeclared_filter_ids(records, {"LowQ"}))

    def test_reported_once_in_first_seen_order(self):
        records = [_record("b"), _record("a;b"), _record("a")]
        self.assertEqual(["b", "a"], scan_undeclared_filter_ids(records, set()))

    def test_capped_at_remaining_vcf_filter_codes(self):
        """ A pathological file mustn't emit thousands of header lines - VCFFilter only has so many
            single-character codes to give out """
        num_codes = VCFFilter.ASCII_MAX - VCFFilter.ASCII_MIN
        records = [_record(f"filter_{i}") for i in range(num_codes + 50)]
        self.assertEqual(num_codes, len(scan_undeclared_filter_ids(records, set())))

    def test_no_codes_left(self):
        declared = {f"declared_{i}" for i in range(VCFFilter.ASCII_MAX - VCFFilter.ASCII_MIN)}
        self.assertEqual([], scan_undeclared_filter_ids([_record("new")], declared))

    def test_short_lines_ignored(self):
        self.assertEqual([], scan_undeclared_filter_ids(["chr1\t100\t.\tT\tA\n"], set()))
