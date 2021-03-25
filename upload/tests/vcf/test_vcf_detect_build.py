from django.conf import settings
from django.test import TestCase
import os
from unittest import skip

from upload.vcf.vcf_import import vcf_detect_genome_build, GenomeBuildDetectionException


class TestVCFDetectBuild(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf", "detect_build_by_header")

    def test_no_contigs(self):
        """ No way to tell what to do """
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "no_contigs.vcf")
        try:
            vcf_detect_genome_build(vcf_filename)
            self.fail("Should have thrown exception for no contigs!")
        except GenomeBuildDetectionException:
            pass

    def test_bad_contigs(self):
        """ We got a VCF from Centogene that had maxint32 for all contig lengths"""
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "bad_contigs.vcf")
        try:
            vcf_detect_genome_build(vcf_filename)
            self.fail("Should have thrown exception for BAD contigs!")
        except GenomeBuildDetectionException:
            pass

    @skip  # TODO: Issue #1857 - Proper handling of hg19 vs GRCh37
    def test_detect_hg19(self):
        """ hg19 - due to MT size  """
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "hg19_contigs.vcf")
        genome_build = vcf_detect_genome_build(vcf_filename)
        self.assertEqual("hg19", genome_build.name, "Matched hg19 genome")

    def test_detect_grch37(self):
        """ GRCh37 - due to MT size  """
        for filename in ["grch37_research_contigs.vcf", "grch37_research_contigs_assembly.vcf"]:
            vcf_filename = os.path.join(self.TEST_DATA_DIR, filename)
            genome_build = vcf_detect_genome_build(vcf_filename)
            self.assertEqual("GRCh37", genome_build.name, "Matched GRCh37 genome")

    def test_freebayes_grch37(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "freebayes_b37.vcf")
        genome_build = vcf_detect_genome_build(vcf_filename)
        self.assertEqual("GRCh37", genome_build.name, "Matched GRCh37 genome")

    def test_big_grch38(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "grch38_huge_header.vcf.gz")
        genome_build = vcf_detect_genome_build(vcf_filename)
        self.assertEqual("GRCh38", genome_build.name, "Matched GRCh38 genome")
