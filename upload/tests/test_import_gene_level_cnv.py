"""Import of a whole-gene CNV VCF - the VCF the loader rewrites, and the events made from it."""
import os

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from genes.gene_copy_number import create_gene_copy_number_events_for_variants
from genes.gene_fusions import create_gene_fusions_for_variants
from genes.models import GeneCopyNumberEvent, GeneCopyNumberEventKind, GeneLevelId
from genes.tests.gene_level_test_utils import create_gene_copy_number_event
from genes.tests.test_gene_fusions import GeneFusionTestCase
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.genomics.vcf_utils import (
    cyvcf2_header_types,
    vcf_get_ref_alt_svlen_and_modification,
)
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from snpdb.models import ImportSource, Variant
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    SimpleVCFImportInfo,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.import_gene_level_cnv_task import (
    GENE_COPY_NUMBER_INFO,
    NO_SEGMENT_MESSAGE,
    GeneLevelCNVCreateVCFTask,
    can_process_file,
    read_gene_level_cnv_records,
)
from upload.vcf.vcf_import import get_copy_number_field, get_gene_level_segment_field

TSO500_DNA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                              "ExampleSample_2600000001", "ExampleSample_DNA_2600000001C")
TSO500_CNV_VCF = os.path.join(TSO500_DNA_DIR, "ExampleSample_DNA_2600000001C.cnv.vcf")
DRAGEN_EXON_CNV_VCF = os.path.join(TSO500_DNA_DIR, "ExampleSample_DNA_2600000001C_DragenExonCNV.vcf")
GENE_LEVEL_CNV_VCF = os.path.join(settings.BASE_DIR, "upload", "test_data", "vcf", "gene_level_cnv",
                                  "gene_level_cnv.vcf")


class TestGeneLevelCNVFactory(TestCase):

    def test_claims_a_segment_field_vcf_over_a_plain_one(self):
        factories = {type(f).__name__: f for f in get_import_task_factories()}
        gene_level = factories["GeneLevelCNVImportTaskFactory"]
        genotype = factories["GenotypeVCFImportFactory"]
        user = User.objects.get_or_create(username='testuser')[0]
        self.assertGreater(gene_level.get_processing_ability(user, TSO500_CNV_VCF, "vcf"),
                           genotype.get_processing_ability(user, TSO500_CNV_VCF, "vcf"))

    def test_a_gene_named_on_exon_level_calls_is_not_a_segment(self):
        """ DragenExonCNV writes GENE=BRCA1 on partial-gene calls - those stay coordinate SVs """
        self.assertFalse(can_process_file(DRAGEN_EXON_CNV_VCF))

    def test_every_called_record_of_the_vendor_file_becomes_an_event(self):
        """ The Illumina example: 16 of its 25 records are called, and none of its gene names are
            ones this test database knows - an unplaceable name still gets an identity """
        reader = cyvcf2.VCF(TSO500_CNV_VCF)
        rewrite = read_gene_level_cnv_records(reader, "SEGID")
        self.assertEqual(25, rewrite.records_read)
        self.assertEqual(16, len(rewrite.records))
        self.assertEqual(0, rewrite.no_segment)
        self.assertEqual(0, rewrite.duplicates)


class TestGeneLevelCNVRewrite(GeneFusionTestCase):
    """ The VCF the loader writes - what the standard insert pipeline then consumes """

    def setUp(self):
        super().setUp()
        user = User.objects.get_or_create(username='testuser')[0]
        file_upload = FileUpload.objects.create(path=GENE_LEVEL_CNV_VCF,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=user,
                                                name="gene_level_cnv.vcf",
                                                file_type=UploadedFileTypes.GENE_LEVEL_CNV_VCF)
        self.upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "gene_level_cnv_variants.vcf")
        self.upload_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                     name="Create Gene-Level CNV VCF", sort_order=0,
                                                     input_filename=GENE_LEVEL_CNV_VCF,
                                                     output_filename=self.vcf_filename)
        self.records_read = GeneLevelCNVCreateVCFTask.process_items(self.upload_step)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        self.records = list(self.reader)
        self.events = {record.INFO.get(GENE_COPY_NUMBER_INFO): record for record in self.records}

    def test_records_are_on_the_gene_level_contig(self):
        self.assertTrue(self.records)
        for record in self.records:
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            _ref, alt, svlen, _mod = vcf_get_ref_alt_svlen_and_modification(
                record, ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self.assertIsNotNone(GeneLevelSymbolicAlt.parse(alt))
            self.assertEqual(0, svlen)

    def test_gain_and_loss_carry_the_gene_in_position_and_alt(self):
        record = self.events["EGFR amplification"]
        egfr = self.hgnc_ids["EGFR"]
        self.assertEqual(egfr, record.POS)
        self.assertEqual(f"<{GeneLevelSymbolicAlt.GAIN}:HGNC:{egfr}>", record.ALT[0])

    def test_no_call_records_are_not_events(self):
        """ BRAF's row is a '.' alt - the caller looked and found nothing """
        self.assertNotIn("BRAF amplification", self.events)
        self.assertNotIn("BRAF loss", self.events)

    def test_a_record_naming_no_gene_is_dropped_and_counted(self):
        self.assertEqual(1, SimpleVCFImportInfo.objects.get(message_string=NO_SEGMENT_MESSAGE).count)
        self.assertEqual(5, self.records_read)
        self.assertEqual(3, len(self.records), "the no-call and the unnamed record are not written")

    def test_an_old_symbol_resolves_to_the_approved_one(self):
        """ The file says SEPT14; the event is SEPTIN14 loss """
        record = self.events["SEPTIN14 loss"]
        self.assertEqual(self.hgnc_ids["SEPTIN14"], record.POS)

    def test_a_name_hgnc_lacks_still_becomes_an_event(self):
        record = self.events["RP11-458D21.5 amplification"]
        gene_level_id = GeneLevelId.objects.get(symbol_str="RP11-458D21.5")
        self.assertGreaterEqual(gene_level_id.pk, GeneLevelId.CUSTOM_ID_START)
        self.assertEqual(gene_level_id.pk, record.POS)

    def test_the_copy_ratio_survives_as_the_copy_number_field(self):
        """ SM is copied through in FORMAT, so the grid's copy number column binds to it """
        header_types = cyvcf2_header_types(self.reader)
        self.assertEqual("SM", get_copy_number_field(set(header_types["FORMAT"]),
                                                     set(header_types.get("INFO", {})),
                                                     single_sample=True))
        self.assertEqual(4.31428, round(self.events["EGFR amplification"].format("SM").flatten()[0], 5))

    def test_the_segment_field_is_declared_so_a_reload_binds_it(self):
        header_types = cyvcf2_header_types(self.reader)
        self.assertEqual("SEGID", get_gene_level_segment_field(set(header_types.get("INFO", {}))))

    def test_the_caller_header_is_kept_so_the_build_is_still_detectable(self):
        self.assertIn("##contig=<ID=chr7,length=159138663>", self.reader.raw_header)


class TestGeneLevelCNVInsert(GeneFusionTestCase):
    """ The post-insert step - GeneCopyNumberEvent rows read off the variants the pipeline inserted """

    def test_creates_an_event_per_gene_level_variant(self):
        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        variant = event.variant
        GeneCopyNumberEvent.objects.filter(pk=event.pk).delete()

        self.assertEqual(1, create_gene_copy_number_events_for_variants(Variant.objects.filter(pk=variant.pk)))
        recreated = GeneCopyNumberEvent.objects.get(variant=variant)
        self.assertEqual(GeneLevelId.objects.get(symbol_str="EGFR"), recreated.gene)
        self.assertEqual("EGFR amplification", recreated.canonical_str)

    def test_is_idempotent(self):
        """ Runs once per pipeline, and a retry runs it again """
        create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        self.assertEqual(0, create_gene_copy_number_events_for_variants(Variant.objects.all()))

    def test_the_fusion_insert_leaves_a_copy_number_variant_alone(self):
        """ Both kinds share the contig, so each claims only its own alts """
        create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        self.assertEqual(0, create_gene_fusions_for_variants(Variant.objects.all()))
