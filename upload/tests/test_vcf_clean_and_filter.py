"""
Tests for the 'vcf_clean_and_filter' management command's sorted check (#127) - the marker file it writes
so preprocess_vcf can tell an unsorted VCF apart from any other failure, and the --allow-unsorted pass
that runs downstream of 'bcftools sort'.

The command only needs the fasta *index* (to rename chroms to fasta sequence names), so these tests point
the build at a placeholder fasta with a .fai built from the build's own contigs.
"""
import copy
import os
import tempfile
from contextlib import redirect_stdout
from io import StringIO

from django.conf import settings
from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase, override_settings

from library.genomics.vcf_utils import write_cleaned_vcf_header
from snpdb.models import GenomeBuild
from upload.management.commands.vcf_clean_and_filter import REFERENCE_SPAN_SKIP_REASON

_TEST_DIR = tempfile.mkdtemp(prefix="test_vcf_clean_and_filter_")
_TEST_FASTA = os.path.join(_TEST_DIR, "placeholder.fna")

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

SORTED_RECORDS = [("1", 200000), ("1", 154560345), ("1", 154560359), ("2", 100000)]
POSITION_BACKWARDS_RECORDS = [("1", 154560359), ("1", 154560345), ("2", 100000)]
UNGROUPED_CONTIG_RECORDS = [("1", 100000), ("2", 100000), ("1", 200000)]


def _write_placeholder_fasta(genome_build):
    with open(_TEST_FASTA, "w") as f:
        f.write(">placeholder\nN\n")  # Sequence is never read - only the index below
    with open(_TEST_FASTA + ".fai", "w") as f:
        for accession, length in genome_build.standard_contigs.values_list("refseq_accession", "length"):
            f.write(f"{accession}\t{length}\t0\t60\t61\n")


def _fake_annotation_settings() -> dict:
    annotation = copy.deepcopy(settings.ANNOTATION)
    annotation[settings.BUILD_GRCH37]["reference_fasta"] = _TEST_FASTA
    return annotation


@override_settings(ANNOTATION=_fake_annotation_settings())
class TestVCFCleanAndFilterSortedCheck(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        _write_placeholder_fasta(cls.genome_build)

    def setUp(self):
        super().setUp()
        self.marker_file = os.path.join(_TEST_DIR, f"{self.id()}.unsorted.txt")

    def _write_vcf(self, records) -> str:
        vcf_filename = os.path.join(_TEST_DIR, f"{self.id()}.vcf")
        with open(vcf_filename, "w") as f:
            f.write(VCF_HEADER)
            for chrom, position in records:
                f.write(f"{chrom}\t{position}\t.\tA\tT\t.\tPASS\t.\n")
        return vcf_filename

    def _run(self, records, **kwargs) -> str:
        """ Returns the written records (header stripped) """
        stdout = StringIO()
        with redirect_stdout(stdout):
            call_command("vcf_clean_and_filter", vcf=self._write_vcf(records),
                         genome_build=self.genome_build.name, **kwargs)
        return "".join(line for line in stdout.getvalue().splitlines(keepends=True) if not line.startswith("#"))

    def test_sorted_vcf_writes_records_and_no_marker(self):
        records = self._run(SORTED_RECORDS, unsorted_marker_file=self.marker_file)
        self.assertEqual(len(records.splitlines()), len(SORTED_RECORDS))
        self.assertFalse(os.path.exists(self.marker_file))

    def test_position_backwards_writes_marker(self):
        with self.assertRaises(CommandError):
            self._run(POSITION_BACKWARDS_RECORDS, unsorted_marker_file=self.marker_file)

        self.assertTrue(os.path.exists(self.marker_file))
        with open(self.marker_file) as f:
            reason = f.read()
        self.assertIn("1:154560345", reason)
        self.assertIn("1:154560359", reason)

    def test_ungrouped_contigs_writes_marker(self):
        with self.assertRaises(CommandError):
            self._run(UNGROUPED_CONTIG_RECORDS, unsorted_marker_file=self.marker_file)

        self.assertTrue(os.path.exists(self.marker_file))
        with open(self.marker_file) as f:
            self.assertIn("grouped together", f.read())

    def test_unsorted_fails_without_marker_file_option(self):
        with self.assertRaises(CommandError):
            self._run(POSITION_BACKWARDS_RECORDS)
        self.assertFalse(os.path.exists(self.marker_file))

    def test_allow_unsorted_writes_all_records(self):
        records = self._run(POSITION_BACKWARDS_RECORDS, allow_unsorted=True)
        self.assertEqual(len(records.splitlines()), len(POSITION_BACKWARDS_RECORDS))
        self.assertFalse(os.path.exists(self.marker_file))

    def test_allow_unsorted_renames_chroms_to_fasta_names(self):
        records = self._run(UNGROUPED_CONTIG_RECORDS, allow_unsorted=True)
        chroms = [line.split("\t")[0] for line in records.splitlines()]
        genome_fasta = self.genome_build.genome_fasta
        expected = [genome_fasta.convert_chrom_to_fasta_sequence(chrom)
                    for chrom, _ in UNGROUPED_CONTIG_RECORDS]
        self.assertEqual(chroms, expected)


@override_settings(ANNOTATION=_fake_annotation_settings())
class TestUndeclaredFiltersSurviveCleanedHeader(TestCase):
    """ preprocess writes the cleaned header then feeds it back as --replace-header, so declaring
        undeclared FILTER values there is what stops vcf_clean_and_filter stripping them """

    FILTER_COLUMN = "LowQ;LowUniqueAlignments"

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        _write_placeholder_fasta(cls.genome_build)

    def _write_source_vcf(self) -> str:
        vcf_filename = os.path.join(_TEST_DIR, f"{self.id()}.vcf")
        with open(vcf_filename, "w") as f:
            f.write(VCF_HEADER.replace(
                "#CHROM", '##FILTER=<ID=LowQ,Description="Below the passing threshold.">\n#CHROM'))
            f.write(f"1\t100000\t.\tT\tA\t0\t{self.FILTER_COLUMN}\t.\n")
        return vcf_filename

    def _cleaned_header(self, vcf_filename) -> str:
        header_filename = os.path.join(_TEST_DIR, f"{self.id()}.header.vcf")
        write_cleaned_vcf_header(self.genome_build, vcf_filename, header_filename)
        with open(header_filename) as f:
            return f.read()

    def test_cleaned_header_declares_undeclared_filter(self):
        header = self._cleaned_header(self._write_source_vcf())
        self.assertIn('##FILTER=<ID=LowUniqueAlignments,Description="Undeclared in source header">', header)
        self.assertIn("##FILTER=<ID=LowQ,", header)

    def test_records_keep_their_filters(self):
        vcf_filename = self._write_source_vcf()
        header_filename = os.path.join(_TEST_DIR, f"{self.id()}.header.vcf")
        write_cleaned_vcf_header(self.genome_build, vcf_filename, header_filename)
        skipped_filters_filename = os.path.join(_TEST_DIR, f"{self.id()}.skipped_filters.txt")

        stdout = StringIO()
        with redirect_stdout(stdout):
            call_command("vcf_clean_and_filter", vcf=vcf_filename, genome_build=self.genome_build.name,
                         replace_header=header_filename, allow_unsorted=True,
                         skipped_filters_stats_file=skipped_filters_filename)

        records = [line for line in stdout.getvalue().splitlines() if not line.startswith("#")]
        self.assertEqual(1, len(records))
        self.assertEqual(self.FILTER_COLUMN, records[0].split("\t")[6])
        self.assertFalse(os.path.exists(skipped_filters_filename))


@override_settings(ANNOTATION=_fake_annotation_settings())
class TestVCFCleanAndFilterReferenceSpans(TestCase):
    """ A reference call over a span carries no call and loses its span on import, so it's skipped -
        the same statement the importer already makes about gVCF reference blocks. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        _write_placeholder_fasta(cls.genome_build)

    def _run(self, records) -> tuple[str, dict]:
        """ Returns (written records with header stripped, skipped record counts) """
        vcf_filename = os.path.join(_TEST_DIR, f"{self.id()}.vcf")
        with open(vcf_filename, "w") as f:
            f.write(VCF_HEADER)
            for alt, info in records:
                f.write(f"1\t100000\t.\tN\t{alt}\t.\tPASS\t{info}\n")

        stats_filename = os.path.join(_TEST_DIR, f"{self.id()}.skipped.txt")
        stdout = StringIO()
        with redirect_stdout(stdout):
            call_command("vcf_clean_and_filter", vcf=vcf_filename, genome_build=self.genome_build.name,
                         allow_unsorted=True, skipped_records_stats_file=stats_filename)

        written = "".join(line for line in stdout.getvalue().splitlines(keepends=True)
                          if not line.startswith("#"))
        skipped = {}
        if os.path.exists(stats_filename):
            with open(stats_filename) as f:
                for line in f:
                    name, count = line.rstrip("\n").split("\t")
                    skipped[name] = int(count)
        return written, skipped

    def test_copy_neutral_segment_skipped(self):
        """ cnv.vcf's copy-neutral rows - ALT='.' over a span with no SVTYPE """
        written, skipped = self._run([(".", "END=104714;REFLEN=4715;SEGID=TNFRSF14")])
        self.assertEqual(written, "")
        self.assertEqual(skipped, {REFERENCE_SPAN_SKIP_REASON: 1})

    def test_undetermined_cnv_call_kept(self):
        """ _DragenExonCNV.vcf's 'Undetermined' BRCA2 row - SVTYPE means the caller is describing an
            event rather than a segmentation interval, so 'we looked and could not tell' survives """
        written, skipped = self._run([(".", "SVTYPE=CNV;END=104714;GENE=BRCA2")])
        self.assertEqual(len(written.splitlines()), 1)
        self.assertEqual(skipped, {})

    def test_single_position_reference_record_kept(self):
        """ Reference variants are first class - only a reference call over a *span* is dropped """
        written, skipped = self._run([(".", ".")])
        self.assertEqual(len(written.splitlines()), 1)
        self.assertEqual(skipped, {})

    def test_symbolic_alts_unaffected(self):
        written, skipped = self._run([("<DUP>", "SVTYPE=CNV;END=104714"),
                                      ("<DEL>", "SVTYPE=CNV;END=104714")])
        self.assertEqual(len(written.splitlines()), 2)
        self.assertEqual(skipped, {})
