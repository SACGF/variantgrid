import tempfile

import cyvcf2
from django.test import TestCase

from library.genomics.vcf_utils import vcf_allele_is_symbolic
from snpdb.variants_to_vcf import VARIANT_GRID_INFO_DICT, _write_sorted_values_to_vcf_file
from snpdb.vcf_export_utils import get_vcf_header_lines


def _sorted_values_from_vcf(filename) -> list[dict]:
    """ Build the `.values()`-shaped dicts the queryset writer consumes, from a real VCF.
        A synthetic 'id' (the record index) stands in for the Variant PK (INFO variant_id). """
    values = []
    for i, v in enumerate(cyvcf2.Reader(filename)):
        svlen = v.INFO.get("SVLEN")
        if isinstance(svlen, tuple):
            svlen = svlen[0]
        end = v.INFO.get("END")
        values.append({
            "locus__contig__name": v.CHROM,
            "locus__position": v.POS,
            "locus__ref__seq": v.REF,
            "alt__seq": v.ALT[0],
            "end": int(end) if end is not None else v.POS + len(v.REF) - 1,
            "svlen": int(svlen) if svlen is not None else None,
            "id": i,
        })
    return values


class TestVariantsToVCF(TestCase):
    """ Tests for the queryset->VCF writer (Writer 1 in the issue #1068 consolidation). Real VCFs
        are read into `.values()`-shaped dicts, written back out, and re-read with a VCF parser. """

    def _write_and_read(self, sorted_values) -> list:
        header_lines = get_vcf_header_lines(info_dict=VARIANT_GRID_INFO_DICT)
        with tempfile.NamedTemporaryFile(mode="wt", suffix=".vcf", delete=True) as temp_file:
            with open(temp_file.name, "wt") as f:
                count = _write_sorted_values_to_vcf_file(header_lines, sorted_values, f,
                                                         info_dict=VARIANT_GRID_INFO_DICT)
            records = list(cyvcf2.Reader(temp_file.name))
        self.assertEqual(len(sorted_values), count)
        return records

    def test_roundtrip(self):
        """ Plain indel + Manta symbolic SVs survive the write/read cycle (incl. variant_id and
            END/SVLEN/SVTYPE for symbolic alts) """
        for filename in ["upload/test_data/vcf/indel_GRCh37.vcf",
                         "upload/test_data/vcf/symbolic_alt/manta.vcf"]:
            with self.subTest(filename=filename):
                sorted_values = _sorted_values_from_vcf(filename)
                records = self._write_and_read(sorted_values)
                for data, v in zip(sorted_values, records):
                    self.assertEqual(data["locus__contig__name"], v.CHROM)
                    self.assertEqual(data["locus__position"], v.POS)
                    self.assertEqual(data["locus__ref__seq"], v.REF)
                    self.assertEqual(data["alt__seq"], v.ALT[0])
                    self.assertEqual(data["id"], v.INFO.get("variant_id"))
                    if vcf_allele_is_symbolic(data["alt__seq"]):
                        self.assertEqual(data["end"], v.INFO.get("END"))
                        self.assertEqual(data["svlen"], v.INFO.get("SVLEN"))
                        self.assertEqual(data["alt__seq"][1:-1], v.INFO.get("SVTYPE"))

    def test_tolerates_missing_end_svlen_keys(self):
        """ Non-symbolic records must not require the 'end'/'svlen' keys (they aren't always
            selected on the alphabetical path). Regression guard for the previous unconditional
            data['end'] lookup that raised KeyError. """
        sorted_values = [
            {"locus__contig__name": "1", "locus__position": 100, "locus__ref__seq": "A",
             "alt__seq": "T", "id": 42},  # no 'end'/'svlen' keys at all
        ]
        records = self._write_and_read(sorted_values)
        self.assertEqual(42, records[0].INFO.get("variant_id"))
