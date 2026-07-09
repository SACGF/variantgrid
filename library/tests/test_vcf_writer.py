import io

from django.test import TestCase

from library.genomics.vcf_writer import (
    VCFInfoHeader,
    build_header_lines,
    symbolic_alt_info,
    VCFWriter,
)

BASE_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


def _rows(text):
    """ Split VCF text into a list of per-line column lists (avoids embedding tabs in tests) """
    return [line.split("\t") for line in text.splitlines()]


class TestVCFWriter(TestCase):
    def test_info_header_str(self):
        h = VCFInfoHeader(id="SVTYPE", type="String", description="Type of structural variant")
        self.assertEqual(
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            str(h))

    def test_info_header_number_and_quote_sanitisation(self):
        h = VCFInfoHeader(id="END", type="Integer", number=".",
                          description='Stop "position" of the interval')
        self.assertEqual(
            "##INFO=<ID=END,Number=.,Type=Integer,Description=\"Stop 'position' of the interval\">",
            str(h))

    def test_build_header_lines_minimal(self):
        lines = build_header_lines()
        self.assertEqual(["##fileformat=VCFv4.1"], lines[:-1])
        self.assertEqual(["#CHROM"] + BASE_COLUMNS[1:], lines[-1].split("\t"))

    def test_build_header_lines_full(self):
        info = [
            VCFInfoHeader(id="END", type="Integer", number=".", description="Stop position of the interval"),
        ]
        formats = ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
        contig_lines = ["##contig=<ID=1,length=249250621>"]
        lines = build_header_lines(meta_lines=["##source=VariantGrid"], info=info, formats=formats,
                                   contig_lines=contig_lines, samples=["S1", "S2"])
        self.assertEqual([
            "##fileformat=VCFv4.1",
            "##source=VariantGrid",
            '##INFO=<ID=END,Number=.,Type=Integer,Description="Stop position of the interval">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=1,length=249250621>",
        ], lines[:-1])
        self.assertEqual(["#CHROM"] + BASE_COLUMNS[1:] + ["FORMAT", "S1", "S2"], lines[-1].split("\t"))

    def test_build_header_lines_samples_without_formats_have_no_format_column(self):
        lines = build_header_lines(samples=["S1"])
        self.assertEqual(["##fileformat=VCFv4.1"], lines[:-1])
        self.assertEqual(["#CHROM"] + BASE_COLUMNS[1:], lines[-1].split("\t"))

    def test_symbolic_alt_info(self):
        self.assertEqual({}, symbolic_alt_info("T", svlen=None, end=100))
        self.assertEqual({"END": 250, "SVLEN": 50, "SVTYPE": "DEL"},
                         symbolic_alt_info("<DEL>", svlen=50, end=250))

    def test_write_record_text_sink(self):
        handle = io.StringIO()
        writer = VCFWriter(handle)
        writer.write_record("1", 100, "A", "T")
        writer.write_record("1", 200, "A", "<DEL>", vcf_id="rs1",
                            info={"SVLEN": 50, "SVTYPE": "DEL", "END": 250})
        self.assertEqual([
            ["1", "100", ".", "A", "T", ".", ".", "."],
            ["1", "200", "rs1", "A", "<DEL>", ".", ".", "SVLEN=50;SVTYPE=DEL;END=250"],
        ], _rows(handle.getvalue()))

    def test_write_to_binary_handle_via_text_wrapper(self):
        # Callers with a binary destination wrap it in io.TextIOWrapper - VCFWriter still only writes str
        raw = io.BytesIO()
        handle = io.TextIOWrapper(raw, encoding="utf-8", write_through=True)
        writer = VCFWriter(handle, header_lines=["##fileformat=VCFv4.1"])
        writer.write_record("1", 100, "A", "T")
        handle.flush()
        lines = raw.getvalue().decode().splitlines()
        self.assertEqual("##fileformat=VCFv4.1", lines[0])
        self.assertEqual(["1", "100", ".", "A", "T", ".", ".", "."], lines[1].split("\t"))

    def test_write_record_with_format_and_samples(self):
        handle = io.StringIO()
        writer = VCFWriter(handle)
        writer.write_record("1", 100, "A", "T", vcf_id="42", info=None,
                            fmt="GT:AD:DP", sample_calls=["0/1:5,6:11", "1/1:0,9:9"])
        self.assertEqual([
            ["1", "100", "42", "A", "T", ".", ".", ".", "GT:AD:DP", "0/1:5,6:11", "1/1:0,9:9"],
        ], _rows(handle.getvalue()))

    def test_encode_info_hook(self):
        handle = io.StringIO()
        writer = VCFWriter(handle, encode_info=lambda v: str(v).replace(";", ",:").replace(",", "|"))
        writer.write_record("1", 100, "A", "T", info={"X": "a;b,c"})
        # ';' -> ',:' then ',' -> '|' (the second replace also hits the ',' just inserted)
        self.assertEqual([["1", "100", ".", "A", "T", ".", ".", "X=a|:b|c"]], _rows(handle.getvalue()))
