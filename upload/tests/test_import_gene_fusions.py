"""Import of AllFusions.csv - the VCF the loader writes, and the GeneFusions made from it."""
import os

import cyvcf2
import simplejson

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.genomics.vcf_utils import vcf_get_ref_alt_svlen_and_modification
from library.genomics.vcf_writer import percent_decode_info_value
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from genes.models import GeneFusion, FusionGeneId
from genes.gene_fusions import create_gene_fusions_for_variants
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.tests.test_gene_fusions import GeneFusionTestCase
from snpdb.models import GenomeBuild, ImportSource, Variant
from upload.gene_fusions.all_fusions_parser import can_process_file, read_all_fusions
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    UploadedFileTypes,
    UploadPipeline,
)
from upload.tasks.import_gene_fusions_task import (
    FUSION_INFO,
    FUSION_OBSERVATIONS_INFO,
    NO_GENOTYPE_CALL,
    GeneFusionCreateVCFTask,
)
from upload.models import UploadStep

TSO500_RNA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                              "ExampleSample_2600000001", "ExampleSample_RNA_2600000001B")
ALL_FUSIONS_CSV = os.path.join(TSO500_RNA_DIR, "ExampleSample_RNA_2600000001B_AllFusions.csv")
EXPECTED_ROWS = 33


class TestAllFusionsParser(TestCase):

    def test_reads_every_row(self):
        comments, rows = read_all_fusions(ALL_FUSIONS_CSV)
        self.assertTrue(comments[0].startswith("# Source = FusionProcessor"))
        self.assertEqual(EXPECTED_ROWS, len(rows))

    def test_ingests_unfiltered(self):
        """ 1 of 149 rows in a real run passes the caller's own filter - we take everything """
        _comments, rows = read_all_fusions(ALL_FUSIONS_CSV)
        filters = {r.data.get("Filter") for r in rows}
        self.assertTrue(any(f and "FAIL" in f for f in filters), "FAIL rows are kept")

    def test_carries_both_callers(self):
        _comments, rows = read_all_fusions(ALL_FUSIONS_CSV)
        self.assertEqual({"DRAGEN", "SpliceGirl"}, {r.caller for r in rows})

    def test_missing_values_are_none(self):
        _comments, rows = read_all_fusions(ALL_FUSIONS_CSV)
        self.assertTrue(any(r.data.get("Score") is None for r in rows), "'N/A' becomes None")

    def test_claims_the_file_over_a_gene_list(self):
        factories = {type(f).__name__: f for f in get_import_task_factories()}
        fusions = factories["GeneFusionsImportTaskFactory"]
        gene_list = factories["GeneListImportTaskFactory"]
        user = User.objects.get_or_create(username='testuser')[0]
        self.assertGreater(fusions.get_processing_ability(user, ALL_FUSIONS_CSV, "csv"),
                           gene_list.get_processing_ability(user, ALL_FUSIONS_CSV, "csv"))

    def test_does_not_claim_other_csvs(self):
        self.assertFalse(can_process_file(os.path.join(TSO500_RNA_DIR,
                                                       "ExampleSample_RNA_2600000001B_SpliceVariants.vcf")))


class TestGeneFusionVCF(GeneFusionTestCase):
    """ The VCF the loader writes - what the standard insert pipeline then consumes """

    def setUp(self):
        super().setUp()
        user = User.objects.get_or_create(username='testuser')[0]
        self.file_upload = FileUpload.objects.create(path=ALL_FUSIONS_CSV,
                                                     import_source=ImportSource.COMMAND_LINE,
                                                     user=user,
                                                     name="ExampleSample_RNA_2600000001B_AllFusions.csv",
                                                     file_type=UploadedFileTypes.GENE_FUSIONS,
                                                     # The file's breakpoints are chrN:pos with nothing
                                                     # to detect a build from, so it has to be declared
                                                     metadata={"genome_build": "GRCh37"})
        self.upload_pipeline = UploadPipeline.objects.create(file_upload=self.file_upload)
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "gene_fusion_variants.vcf")
        create_vcf_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                    name="Create Gene Fusion Variant VCF", sort_order=0,
                                                    input_filename=ALL_FUSIONS_CSV,
                                                    output_filename=self.vcf_filename)
        self.rows_processed = GeneFusionCreateVCFTask.process_items(create_vcf_step)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        self.records = list(self.reader)

    def test_all_rows_read(self):
        self.assertEqual(EXPECTED_ROWS, self.rows_processed)

    def test_header_declares_the_sample_and_source(self):
        """ What ImportCreateVCFModelForGenotypeVCFTask makes the VCF and Sample from """
        self.assertEqual([self.file_upload.name], self.reader.samples)
        self.assertIn("FusionProcessor 1.0.0.614", self.reader.raw_header)

    def test_records_are_on_the_gene_level_contig(self):
        self.assertTrue(self.records)
        for record in self.records:
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            _ref, alt, svlen, _mod = vcf_get_ref_alt_svlen_and_modification(
                record, ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self.assertIsNotNone(GeneLevelSymbolicAlt.parse(alt))
            self.assertEqual(0, svlen)

    def test_genotype_asserts_presence_not_a_diploid_call(self):
        for record in self.records:
            self.assertIn(NO_GENOTYPE_CALL, str(record))

    def test_sept14_resolves_to_septin14(self):
        """ The file says SEPT14; the fusion is EGFR::SEPTIN14 """
        fusions = {record.INFO.get(FUSION_INFO) for record in self.records}
        self.assertIn("EGFR::SEPTIN14", fusions)

    def test_repeated_gene_pair_collapses_keeping_every_observation(self):
        """ ENTPD3::RPL14 appears three times from one caller with three different 5' breakpoints """
        by_fusion = {record.INFO.get(FUSION_INFO): record for record in self.records}
        record = by_fusion["ENTPD3::RPL14"]
        observations = simplejson.loads(percent_decode_info_value(record.INFO.get(FUSION_OBSERVATIONS_INFO)))
        self.assertEqual(3, len(observations))
        self.assertEqual(3, len({o["Gene A Breakpoint"] for o in observations}))


class TestGeneFusionInsert(GeneFusionTestCase):
    """ The post-insert step - GeneFusion rows read off the variants the pipeline inserted """

    def test_creates_a_gene_fusion_per_gene_level_variant(self):
        gene_fusion = create_gene_fusion("BCR", "ABL1")
        variant = gene_fusion.variant
        GeneFusion.objects.filter(pk=gene_fusion.pk).delete()

        self.assertEqual(1, create_gene_fusions_for_variants(Variant.objects.filter(pk=variant.pk)))
        recreated = GeneFusion.objects.get(variant=variant)
        self.assertEqual(FusionGeneId.objects.get(symbol_str="BCR"), recreated.anchor)
        self.assertEqual(FusionGeneId.objects.get(symbol_str="ABL1"), recreated.partner)
        self.assertTrue(recreated.is_ordered)

    def test_is_idempotent(self):
        """ Runs once per pipeline, and a retry runs it again """
        create_gene_fusion("BCR", "ABL1")
        self.assertEqual(0, create_gene_fusions_for_variants(Variant.objects.all()))
