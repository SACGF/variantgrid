import os

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.fake_annotation import get_fake_annotation_version
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import sha256sum_str
from snpdb.models import ImportSource, Sequence
from upload.models import (
    FileUpload,
    ModifiedImportedVariantOperation,
    ModifiedImportedVariants,
    UploadedFileTypes,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
)
from upload.vcf.bulk_genotype_vcf_processor import BulkGenotypeVCFProcessor
from upload.vcf.bulk_no_genotype_vcf_processor import BulkNoGenotypeVCFProcessor
from upload.vcf.sql_copy_files import COHORT_GENOTYPE_HEADER
from upload.vcf.vcf_import import create_cohort_genotype_collection_from_vcf, create_vcf_from_vcf


class TestVCFProcessors(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf")
    TSO500_RNA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                                  "ExampleSample_2600000001", "ExampleSample_RNA_2600000001B")
    TSO500_DNA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                                  "ExampleSample_2600000001", "ExampleSample_DNA_2600000001C")

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_sha256_hash=sha256sum_str(base))

    @classmethod
    def _create_fake_upload_step(cls, vcf_filename, metadata=None) -> UploadStep:
        user = User.objects.get_or_create(username='testuser')[0]
        file_upload = FileUpload.objects.create(path=vcf_filename,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=user,
                                                file_type=UploadedFileTypes.VCF,
                                                metadata=metadata or {})

        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        return UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                         input_filename=vcf_filename,
                                         sort_order=0)

    @classmethod
    def _create_fake_upload_step_and_vcf(cls, vcf_filename, vcf_reader) -> tuple[UploadStep, UploadedVCF]:
        upload_step = cls._create_fake_upload_step(vcf_filename)
        vcf = create_vcf_from_vcf(upload_step, vcf_reader)
        get_fake_annotation_version(vcf.genome_build)
        create_cohort_genotype_collection_from_vcf(vcf, vcf_reader)
        uploaded_vcf = UploadedVCF.objects.get(upload_pipeline=upload_step.upload_pipeline)
        return upload_step, uploaded_vcf

    def _test_genotype_processor(self, vcf_filename, processor_klass):
        """ I keep forgetting to adjust the columns to match the CSV """
        fast_vcf_reader = cyvcf2.VCF(vcf_filename)
        upload_step, uploaded_vcf = self._create_fake_upload_step_and_vcf(vcf_filename, fast_vcf_reader)
        cohort_genotype_collection = uploaded_vcf.vcf.cohort.cohort_genotype_collection
        processor = processor_klass(upload_step, cohort_genotype_collection, uploaded_vcf, None)

        for v in fast_vcf_reader:
            processor.process_entry(v)
            break

        cg = None
        for field_name in ["locus_cohort_genotypes", "cohort_genotypes"]:
            if f := getattr(processor, field_name, None):
                cg = f[0]
                break
        if cg is None:
            raise ValueError("Couldn't find array to retrieve cohort genotype")
        len_genotype_cols = len(cg)
        len_columns = len(COHORT_GENOTYPE_HEADER) - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS
        message = f"{processor_klass} CohortGenotypeData ({len_genotype_cols} cols) !=  CSV columns ({len_columns})"
        self.assertEqual(len_genotype_cols, len_columns, message)

    def test_no_genotype_processor(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "no_genotype.GRCh37.vcf")
        self._test_genotype_processor(vcf_filename, BulkNoGenotypeVCFProcessor)

    def test_genotype_processor(self):
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "sample1_hg19.vcf")
        self._test_genotype_processor(vcf_filename, BulkGenotypeVCFProcessor)

    def _create_vcf_with_metadata(self, vcf_filename, metadata):
        upload_step = self._create_fake_upload_step(vcf_filename, metadata)
        return create_vcf_from_vcf(upload_step, cyvcf2.VCF(vcf_filename))

    def test_declared_genome_build_reaches_the_vcf(self):
        """ _DragenExonCNV.vcf has no contigs and an unresolvable '##reference', so without a declared
            build create_vcf_from_vcf leaves it null and the pipeline stops at REQUIRES_USER_INPUT """
        vcf_filename = os.path.join(self.TSO500_DNA_DIR,
                                    "ExampleSample_DNA_2600000001C_DragenExonCNV.vcf")
        vcf = self._create_vcf_with_metadata(vcf_filename, {"genome_build": "GRCh37"})
        self.assertEqual("GRCh37", vcf.genome_build.name)

    def test_declared_source_reaches_the_vcf(self):
        """ cnv.vcf declares no '##source' at all, so nothing keyed on the caller can reach it """
        vcf_filename = os.path.join(self.TSO500_DNA_DIR, "ExampleSample_DNA_2600000001C.cnv.vcf")
        self.assertEqual("", self._create_vcf_with_metadata(vcf_filename, {}).source)
        vcf = self._create_vcf_with_metadata(vcf_filename, {"source": "DRAGEN TSO500 CNV"})
        self.assertEqual("DRAGEN TSO500 CNV", vcf.source)

    def test_declared_source_wins_over_header(self):
        vcf_filename = os.path.join(self.TSO500_RNA_DIR, "ExampleSample_RNA_2600000001B_SpliceVariants.vcf")
        vcf = self._create_vcf_with_metadata(vcf_filename, {"source": "SpliceGirl"})
        self.assertEqual("SpliceGirl", vcf.source)
        self.assertIn("SpliceGirl 1.0.0.614", vcf.header, "Header text is kept either way")

    def test_splicegirl_sample_field_overrides(self):
        """ SpliceGirl reuses AD for splice-supporting reads and DP for *reference* reads, so binding
            them by name gives empty allele depth, no VAF and a read depth column showing ref counts.
            The ^SpliceGirl VCFSourceSettings row remaps them, making VAF = AD / (AD + DP) - which in
            this caller's files is exactly INFO/ALTDEDUP / (ALTDEDUP + REFDEDUP). """

        vcf_filename = os.path.join(self.TSO500_RNA_DIR, "ExampleSample_RNA_2600000001B_SpliceVariants.vcf")
        fast_vcf_reader = cyvcf2.VCF(vcf_filename)
        _upload_step, uploaded_vcf = self._create_fake_upload_step_and_vcf(vcf_filename, fast_vcf_reader)
        vcf = uploaded_vcf.vcf

        self.assertEqual("SpliceGirl 1.0.0.614", vcf.source)
        self.assertIsNone(vcf.allele_depth_field, "AD isn't the packed ref,alt array cyvcf2 expects")
        self.assertEqual("AD", vcf.alt_depth_field)
        self.assertEqual("DP", vcf.ref_depth_field)
        self.assertIsNone(vcf.read_depth_field, "DP here is reference reads, not total depth")
        self.assertIsNone(vcf.allele_frequency_field, "VAF is derived, not read")

    @override_settings(VARIANT_SYMBOLIC_ALT_SIZE=1)
    def _process_splice_vcf(self):
        """ Runs every record through the processor, returning its (chrom, pos) list alongside it.

            VARIANT_SYMBOLIC_ALT_SIZE keeps the <DEL>s symbolic - a real deployment expands the shorter
            ones against the reference fasta, which is a different question from the VAF being tested """
        for symbolic_alt in [VCFSymbolicAllele.DEL]:
            Sequence.objects.get_or_create(seq=symbolic_alt, seq_sha256_hash=sha256sum_str(symbolic_alt))

        vcf_filename = os.path.join(self.TSO500_RNA_DIR, "ExampleSample_RNA_2600000001B_SpliceVariants.vcf")
        fast_vcf_reader = cyvcf2.VCF(vcf_filename)
        upload_step, uploaded_vcf = self._create_fake_upload_step_and_vcf(vcf_filename, fast_vcf_reader)
        cohort_genotype_collection = uploaded_vcf.vcf.cohort.cohort_genotype_collection
        preprocess_vcf_import_info = ModifiedImportedVariants.objects.create(upload_step=upload_step)
        processor = BulkGenotypeVCFProcessor(upload_step, cohort_genotype_collection, uploaded_vcf,
                                             preprocess_vcf_import_info)

        positions = []
        for v in fast_vcf_reader:
            processor.process_entry(v)
            positions.append((v.CHROM, v.POS))
        processor.finished_locus()
        return positions, processor

    @staticmethod
    def _vaf(cohort_genotype, index) -> float:
        return float(cohort_genotype[index].strip("{}"))

    def test_splicegirl_derived_vaf(self):
        positions, processor = self._process_splice_vcf()
        self.assertEqual(17, len(positions), "All 17 splice records processed")
        self.assertEqual(17, len(processor.cohort_genotypes))

        vaf_index = processor.cohort_gt_vaf_index
        vaf_by_position = {}
        for position, cohort_genotype in zip(positions, processor.cohort_genotypes):
            vaf_by_position.setdefault(position, []).append(self._vaf(cohort_genotype, vaf_index))

        # EGFRvIII: ALTDEDUP=64, REFDEDUP=1
        self.assertAlmostEqual(64 / 65, vaf_by_position[("chr7", 55087058)][0], places=3)
        # MET exon 14 skipping: ALTDEDUP=91, REFDEDUP=1
        self.assertAlmostEqual(91 / 92, vaf_by_position[("chr7", 116411708)][0], places=3)
        # Background call: ALTDEDUP=1, REFDEDUP=80
        self.assertAlmostEqual(1 / 81, vaf_by_position[("chr1", 120464432)][0], places=3)

    def test_splicegirl_allele_depth_is_alt_reads(self):
        positions, processor = self._process_splice_vcf()
        allele_depth_index = COHORT_GENOTYPE_HEADER.index("samples_allele_depth") \
            - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS
        read_depth_index = COHORT_GENOTYPE_HEADER.index("samples_read_depth") \
            - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS

        by_position = dict(zip(positions, processor.cohort_genotypes))
        egfr = by_position[("chr7", 55087058)]
        self.assertEqual("{64}", egfr[allele_depth_index], "Allele depth is AD (alt reads)")
        self.assertIsNone(egfr[read_depth_index], "Read depth not populated from DP (reference reads)")

    def test_splicegirl_shared_locus_is_recorded(self):
        """ chr2:47637511 appears twice with different END. process_entry keys the locus on
            (CHROM, POS, ref), so both join one locus and their depths sum - each record's VAF comes
            out 1/182 rather than 1/91. Accepted, but recorded so it stays reconstructable. """

        positions, processor = self._process_splice_vcf()
        vaf_index = processor.cohort_gt_vaf_index
        shared = [cg for position, cg in zip(positions, processor.cohort_genotypes)
                  if position == ("chr2", 47637511)]
        self.assertEqual(2, len(shared))
        for cohort_genotype in shared:
            self.assertAlmostEqual(1 / 182, self._vaf(cohort_genotype, vaf_index), places=4)

        shared_locus_mivs = [miv for miv in processor.modified_imported_variants
                             if miv[0] == ModifiedImportedVariantOperation.SHARED_LOCUS]
        self.assertEqual(2, len(shared_locus_mivs), "One per record sharing the locus")
        for miv in shared_locus_mivs:
            detail = miv[-1]
            self.assertIn("summed across 2 VCF records", detail)
            self.assertIn("[182]", detail, "Denominator recorded")
            self.assertIn("ref=[90], alt=[1]", detail, "This record's own depths recorded")

    def test_no_shared_locus_records_for_ordinary_vcf(self):
        """ Decomposed multi-allelics are what the summing is *for* and already get a NORMALIZATION
            row, so an ordinary one-record-per-locus VCF records nothing """
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "sample1_hg19.vcf")
        fast_vcf_reader = cyvcf2.VCF(vcf_filename)
        upload_step, uploaded_vcf = self._create_fake_upload_step_and_vcf(vcf_filename, fast_vcf_reader)
        processor = BulkGenotypeVCFProcessor(upload_step, uploaded_vcf.vcf.cohort.cohort_genotype_collection,
                                             uploaded_vcf,
                                             ModifiedImportedVariants.objects.create(upload_step=upload_step))
        for v in fast_vcf_reader:
            processor.process_entry(v)
        processor.finished_locus()

        operations = {miv[0] for miv in processor.modified_imported_variants}
        self.assertNotIn(ModifiedImportedVariantOperation.SHARED_LOCUS, operations)

    def test_genotype_processor_undeclared_filter(self):
        """VCF with a FILTER value not declared in the header should not crash the processor.
        The processor must auto-create a VCFFilter record for the undeclared code.
        """
        vcf_filename = os.path.join(self.TEST_DATA_DIR, "undeclared_filter.vcf")
        self._test_genotype_processor(vcf_filename, BulkGenotypeVCFProcessor)

