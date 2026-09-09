"""Import of AllFusions.csv - the VCF the loader writes, and the GeneFusions made from it."""
import os
import tempfile

import cyvcf2
import simplejson
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from genes.gene_fusions import GeneFusionResolver, create_gene_fusions_for_variants
from genes.models import FusionGeneId, GeneFusion
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.tests.test_gene_fusions import GeneFusionTestCase
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.genomics.vcf_utils import vcf_get_ref_alt_svlen_and_modification
from library.genomics.vcf_writer import percent_decode_info_value
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from snpdb.models import VCF, GenomeBuild, ImportSource, Variant, VCFSourceSettings
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.import_dragen_tso500_all_fusions_task import (
    ALT_READS_FORMAT,
    FUSION_INFO,
    FUSION_OBSERVATIONS_INFO,
    REF_READS_FORMAT,
    DragenTSO500AllFusionsCreateVCFTask,
)
from upload.tso500.dragen_all_fusions_parser import can_process_file, read_all_fusions
from upload.vcf.vcf_import import resolve_genome_build

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
        fusions = factories["DragenTSO500AllFusionsImportTaskFactory"]
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
                                                     file_type=UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS,
                                                     metadata={"genome_build": "GRCh37"})
        self.upload_pipeline = UploadPipeline.objects.create(file_upload=self.file_upload)
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "gene_fusion_variants.vcf")
        create_vcf_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                    name="Create Gene Fusion Variant VCF", sort_order=0,
                                                    input_filename=ALL_FUSIONS_CSV,
                                                    output_filename=self.vcf_filename)
        self.rows_processed = DragenTSO500AllFusionsCreateVCFTask.process_items(create_vcf_step)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        self.records = list(self.reader)

    def test_all_rows_read(self):
        self.assertEqual(EXPECTED_ROWS, self.rows_processed)

    def test_header_declares_the_sample_and_source(self):
        """ What ImportCreateVCFModelForGenotypeVCFTask makes the VCF and Sample from """
        self.assertEqual([self.file_upload.name], self.reader.samples)
        self.assertIn("FusionProcessor 1.0.0.614", self.reader.raw_header)

    def test_genome_build_comes_from_the_source(self):
        """ Neither the csv nor the VCF written from it names a build, and the gene-level contig
            matches every one - the ^FusionProcessor VCFSourceSettings row says which build TSO500
            is run against, so the import doesn't stop at REQUIRES_USER_INPUT """
        self.file_upload.metadata = {}
        self.assertEqual(GenomeBuild.grch37(), resolve_genome_build(self.reader, self.file_upload))

    def test_records_are_on_the_gene_level_contig(self):
        self.assertTrue(self.records)
        for record in self.records:
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            _ref, alt, svlen, _mod = vcf_get_ref_alt_svlen_and_modification(
                record, ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self.assertIsNotNone(GeneLevelSymbolicAlt.parse(alt))
            self.assertEqual(0, svlen)

    def test_sample_carries_read_support_not_a_genotype(self):
        """ A caller asserts the fusion is present, so there is no GT to filter zygosity on - the
            sample column holds how many reads support it and how many don't """
        for record in self.records:
            self.assertEqual([ALT_READS_FORMAT, REF_READS_FORMAT], record.FORMAT)

    def test_read_support_sums_breakpoints_and_keeps_shared_reference(self):
        """ ENTPD3::RPL14's three calls have 4, 12 and 10 supporting reads at distinct 5' breakpoints
            and all re-report the same 1 + 3108 reference reads across the junctions """
        by_fusion = {record.INFO.get(FUSION_INFO): record for record in self.records}
        record = by_fusion["ENTPD3::RPL14"]
        self.assertEqual(26, record.format(ALT_READS_FORMAT).flatten()[0])
        self.assertEqual(3109, record.format(REF_READS_FORMAT).flatten()[0])

    def test_source_settings_bind_read_support_as_depth(self):
        """ The ^FusionProcessor row makes ALT_READS/REF_READS the depths the sample node and VAF
            use, and with no GT in the header nothing binds a genotype """
        vcf = VCF(source="FusionProcessor 1.0.0.614", genotype_samples=1)
        for vss in VCFSourceSettings.get_for_source(vcf.source):
            vss.apply_sample_field_overrides(vcf)
        self.assertEqual(ALT_READS_FORMAT, vcf.alt_depth_field)
        self.assertEqual(REF_READS_FORMAT, vcf.ref_depth_field)
        self.assertTrue(vcf.has_depth)
        self.assertFalse(vcf.has_genotype)

    def test_create_vcf_step_resolves_the_build_from_the_source_line(self):
        """ The step needs a build before the VCF exists, to look breakpoints up in - with nothing
            declared at upload the '# Source =' line is what answers it """
        self.file_upload.metadata = {}
        self.file_upload.save()
        vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "gene_fusion_variants_no_metadata.vcf")
        upload_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                name="Create Gene Fusion Variant VCF", sort_order=1,
                                                input_filename=ALL_FUSIONS_CSV,
                                                output_filename=vcf_filename)
        self.assertEqual(EXPECTED_ROWS, DragenTSO500AllFusionsCreateVCFTask.process_items(upload_step))
        fusions = {record.INFO.get(FUSION_INFO) for record in cyvcf2.VCF(vcf_filename)}
        self.assertIn("EGFR::SEPTIN14", fusions, "the approved symbol, whichever way the build came")

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

    def test_variant_inserted_with_zero_svlen(self):
        """ The pipeline inserts through VariantPKLookup - gene-level variants must land with
            svlen=0, as unique_together does nothing on null @see snpdb.gene_level_variants """
        resolver = GeneFusionResolver()
        resolved_fusion = resolver.resolve_fusion(resolver.resolve_side("BCR"),
                                                  resolver.resolve_side("ABL1"), True)
        variant_coordinate = resolved_fusion.variant_coordinate
        genome_build = GenomeBuild.grch37()
        with tempfile.TemporaryDirectory() as working_dir:
            variant_pk_lookup = VariantPKLookup(genome_build, working_dir=working_dir)
            variant_pk_lookup.add(variant_coordinate)
            variant_pk_lookup.batch_check(insert_unknown=True)

        variant = Variant.get_from_variant_coordinate(variant_coordinate, genome_build)
        self.assertEqual(0, variant.svlen)

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
