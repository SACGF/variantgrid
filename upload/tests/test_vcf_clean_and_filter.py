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

from snpdb.models import GenomeBuild

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
        with open(_TEST_FASTA, "w") as f:
            f.write(">placeholder\nN\n")  # Sequence is never read - only the index below
        with open(_TEST_FASTA + ".fai", "w") as f:
            for accession, length in cls.genome_build.standard_contigs.values_list("refseq_accession", "length"):
                f.write(f"{accession}\t{length}\t0\t60\t61\n")

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
