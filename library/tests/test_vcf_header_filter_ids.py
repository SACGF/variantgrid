from django.test import TestCase

from library.genomics.vcf_utils import vcf_header_filter_ids


class TestVCFHeaderFilterIds(TestCase):
    def test_reads_ids(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##FILTER=<ID=LowQ,Description="Low quality">\n',
            '##FILTER=<ID=base_quality,Description="Site filtered because median base quality is low">\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        self.assertEqual({"LowQ", "base_quality"}, vcf_header_filter_ids(header))

    def test_extra_keys(self):
        """ Illumina's LrCalculator writes Number/Type on FILTER lines - strict parsers reject the whole header """
        header = ['##FILTER=<ID=Undetermined,Number=1,Type=String,Description="CNV type may be undetermined. ">\n']
        self.assertEqual({"Undetermined"}, vcf_header_filter_ids(header))

    def test_id_only(self):
        self.assertEqual({"PASS"}, vcf_header_filter_ids(['##FILTER=<ID=PASS>\n']))

    def test_ignores_other_header_types(self):
        header = [
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            '##ALT=<ID=DEL,Description="Deletion relative to the reference.">\n',
        ]
        self.assertEqual(set(), vcf_header_filter_ids(header))

    def test_no_filters(self):
        self.assertEqual(set(), vcf_header_filter_ids([]))
