import os
import tempfile

from django.test import TestCase

from library.genomics.vcf_utils import get_variant_caller_and_version_from_vcf

COLUMN_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class TestVCFCallerAndVersion(TestCase):
    def _caller_and_version(self, *meta_lines):
        with tempfile.TemporaryDirectory() as tmp_dir:
            filename = os.path.join(tmp_dir, "test.vcf")
            with open(filename, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.writelines(line + "\n" for line in meta_lines)
                f.write(COLUMN_HEADER)
            return get_variant_caller_and_version_from_vcf(filename)

    def test_source_space_and_underscore(self):
        self.assertEqual(("freeBayes", "1.3.5"), self._caller_and_version("##source=freeBayes v1.3.5"))
        self.assertEqual(("VarDict", "1.8.2"), self._caller_and_version("##source=VarDict_v1.8.2"))

    def test_gatk_commandline(self):
        """ HaplotypeCaller is reported as plain GATK, Version quotes are stripped """
        meta = '##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="x --y",Version="4.1.2.0",Date="d">'
        self.assertEqual(("GATK", "4.1.2.0"), self._caller_and_version(meta))

    def test_unrecognised_source(self):
        self.assertEqual((None, None), self._caller_and_version("##source=ClinVar"))
