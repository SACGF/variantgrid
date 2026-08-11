import os
from unittest import skip

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import ImportSource
from upload.models import FileUpload, UploadedFileTypes
from upload.vcf.vcf_import import (
    GenomeBuildDetectionException,
    GenomeBuildMismatchException,
    resolve_genome_build,
    vcf_detect_genome_build,
)


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


class TestResolveGenomeBuild(TestCase):
    """ resolve_genome_build combines header detection with the build declared at upload
        @see upload.upload_metadata """

    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf", "detect_build_by_header")
    TSO500_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                              "ExampleSample_2600000001", "ExampleSample_DNA_2600000001C")

    @staticmethod
    def _file_upload(user, metadata: dict) -> FileUpload:
        return FileUpload.objects.create(user=user, name="declared.vcf", path="declared.vcf",
                                         import_source=ImportSource.COMMAND_LINE,
                                         file_type=UploadedFileTypes.VCF, metadata=metadata)

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]

    def _resolve(self, vcf_filename, metadata=None):
        file_upload = self._file_upload(self.user, metadata) if metadata is not None else None
        return resolve_genome_build(cyvcf2.VCF(vcf_filename), file_upload)

    def test_declared_build_used_when_header_has_none(self):
        """ _DragenExonCNV.vcf has no contigs and an unresolvable '##reference' - it only loads
            because the submitter declared the build """
        vcf_filename = os.path.join(self.TSO500_DIR, "ExampleSample_DNA_2600000001C_DragenExonCNV.vcf")
        self.assertIsNone(self._resolve(vcf_filename, {}))
        genome_build = self._resolve(vcf_filename, {"genome_build": "GRCh37"})
        self.assertEqual("GRCh37", genome_build.name)

    def test_header_detection_still_wins_when_they_agree(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "grch37_research_contigs.vcf")
        genome_build = self._resolve(vcf_filename, {"genome_build": "GRCh37"})
        self.assertEqual("GRCh37", genome_build.name)

    def test_disagreement_fails_rather_than_picking_a_winner(self):
        """ Differing contigs mean the coordinates aren't what the submitter thinks they are """
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "grch37_research_contigs.vcf")
        with self.assertRaises(GenomeBuildMismatchException):
            self._resolve(vcf_filename, {"genome_build": "GRCh38"})
