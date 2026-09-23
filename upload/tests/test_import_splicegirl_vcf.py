"""SpliceGirl's SpliceVariants.vcf rewritten as gene-level splice variants (#1903)."""
import os

import cyvcf2
import simplejson
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from genes.models import HGNC, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import get_sequence, make_release_gene
from library.genomics.vcf_writer import percent_decode_info_value
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from snpdb.models import GenomeBuild, ImportSource
from upload.import_task_factories.import_task_factory import get_import_task_factory_from_extension
from upload.models import (
    FileUpload,
    ModifiedImportedVariants,
    SimpleVCFImportInfo,
    UploadedFileTypes,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.import_splicegirl_vcf_task import (
    UNRESOLVED_MESSAGE,
    SpliceGirlCreateVCFTask,
)
from upload.tso500.dragen_combined_variant_output_parser import (
    SPLICE_INFO,
    SPLICE_OBSERVATION_INFO,
    format_splice_observation,
)
from upload.vcf.bulk_genotype_vcf_processor import BulkGenotypeVCFProcessor
from upload.vcf.sql_copy_files import COHORT_GENOTYPE_HEADER
from upload.vcf.vcf_import import create_cohort_genotype_collection_from_vcf, create_vcf_from_vcf

TSO500_PAIR_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500", "ExampleSample_2600000001")
SPLICE_VARIANTS_VCF = os.path.join(TSO500_PAIR_DIR, "ExampleSample_RNA_2600000001B",
                                   "ExampleSample_RNA_2600000001B_SpliceVariants.vcf")
CNV_VCF = os.path.join(TSO500_PAIR_DIR, "ExampleSample_DNA_2600000001C", "ExampleSample_DNA_2600000001C.cnv.vcf")
RECORDS = 18
SEEDED_RECORDS = 3
# The background junctions this test's release puts in a gene: three in NOTCH2, both chr2:47637511 in
# MSH2 - the other ten are in none
IN_A_GENE_RECORDS = 5


class TestSpliceGirlVCF(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(644, "AR"), (3236, "EGFR"), (7029, "MET"), (7325, "MSH2"), (7882, "NOTCH2")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk
        annotation_version = get_fake_annotation_version(cls.genome_build)
        release = annotation_version.variant_annotation_version.gene_annotation_release
        make_release_gene(cls.genome_build, release, "ENSG00000134250", "NOTCH2",
                          "ENST00000256646.1", "1", 120_460_000, hgnc_id=7882)
        make_release_gene(cls.genome_build, release, "ENSG00000095002", "MSH2",
                          "ENST00000233146.1", "2", 47_630_000, hgnc_id=7325)

    def setUp(self):
        self.user = User.objects.get_or_create(username='testuser')[0]
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "gene_level_splice.vcf")
        self.upload_step = self._create_upload_step()
        self.records_read = SpliceGirlCreateVCFTask.process_items(self.upload_step)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        self.records = list(self.reader)

    def _create_upload_step(self) -> UploadStep:
        file_upload = FileUpload.objects.create(path=SPLICE_VARIANTS_VCF,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=self.user,
                                                name=os.path.basename(SPLICE_VARIANTS_VCF),
                                                file_type=UploadedFileTypes.GENE_LEVEL_SPLICE_VCF)
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        return UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                         name="Create Gene-Level Splice VCF", sort_order=0,
                                         input_filename=SPLICE_VARIANTS_VCF,
                                         output_filename=self.vcf_filename)

    def _by_splice(self) -> dict:
        return {percent_decode_info_value(record.INFO.get(SPLICE_INFO)): record
                for record in self.records}

    def test_upload_picks_the_splice_loader_off_the_header(self):
        """ The pipeline and the upload page both send it as a plain VCF - '##source=SpliceGirl' is
            what routes it, whatever the file is called """
        factory = get_import_task_factory_from_extension(self.user, SPLICE_VARIANTS_VCF, "vcf")
        self.assertEqual(UploadedFileTypes.GENE_LEVEL_SPLICE_VCF, factory.get_uploaded_file_type())
        factory = get_import_task_factory_from_extension(self.user, CNV_VCF, "vcf")
        self.assertNotEqual(UploadedFileTypes.GENE_LEVEL_SPLICE_VCF, factory.get_uploaded_file_type())

    def test_every_junction_in_a_gene_is_kept_with_its_filter(self):
        """ PASS and LowQ alike - the analysis filters decide, not the loader """
        self.assertEqual(RECORDS, self.records_read)
        self.assertEqual(SEEDED_RECORDS + IN_A_GENE_RECORDS, len(self.records))
        self.assertTrue(all(record.CHROM == GENE_LEVEL_CONTIG_NAME for record in self.records))
        filters = [record.FILTER for record in self.records]
        self.assertEqual(SEEDED_RECORDS, filters.count(None), "PASS")
        self.assertEqual(IN_A_GENE_RECORDS, filters.count("LowQ;LowUniqueAlignments"))

    def test_junctions_in_no_gene_are_a_message(self):
        message = SimpleVCFImportInfo.objects.get(upload_step=self.upload_step, message_string=UNRESOLVED_MESSAGE)
        self.assertEqual(RECORDS - SEEDED_RECORDS - IN_A_GENE_RECORDS, message.count)

    def test_seeded_junctions_get_their_label(self):
        """ The same alts the CombinedVariantOutput's rows resolved to, so earlier imports and
            classifications land on one Variant """
        alts = {record.ALT[0] for record in self.records}
        self.assertTrue({f"<SPLICE:HGNC:{self.hgnc_ids['AR']}:V_7>",
                         f"<SPLICE:HGNC:{self.hgnc_ids['EGFR']}:V_III>",
                         f"<SPLICE:HGNC:{self.hgnc_ids['MET']}:EXON_14_SKIPPING>"}.issubset(alts))

    def test_one_donor_two_acceptors_is_two_variants(self):
        """ chr2:47637511 is written twice with different END - as <DEL>s both joined one Locus and
            summed their depths; as junctions they are two events """
        msh2_alts = {record.ALT[0] for record in self.records if record.POS == self.hgnc_ids["MSH2"]}
        self.assertEqual(2, len(msh2_alts))

    def test_the_callers_record_rides_along_in_info(self):
        encoded = self._by_splice()["MET exon 14 skipping"].INFO.get(SPLICE_OBSERVATION_INFO)
        observation = simplejson.loads(percent_decode_info_value(encoded))
        self.assertEqual("chr7:116411708→chr7:116414934 (91 reads)", format_splice_observation(observation))
        self.assertEqual(91, observation["ALTDUP"])

    def _process(self) -> tuple[list, BulkGenotypeVCFProcessor]:
        upload_step = UploadStep.objects.create(upload_pipeline=self.upload_step.upload_pipeline,
                                                input_filename=self.vcf_filename, sort_order=1)
        vcf = create_vcf_from_vcf(upload_step, self.reader)
        create_cohort_genotype_collection_from_vcf(vcf, self.reader)
        records = list(cyvcf2.VCF(self.vcf_filename))
        for record in records:  # The processor reads the Sequence table once, when it is made
            get_sequence(record.REF)
            get_sequence(record.ALT[0])
        uploaded_vcf = upload_step.upload_pipeline.uploadedvcf
        processor = BulkGenotypeVCFProcessor(upload_step, vcf.cohort.cohort_genotype_collection, uploaded_vcf,
                                             ModifiedImportedVariants.objects.create(upload_step=upload_step))
        for record in records:
            processor.process_entry(record)
        processor.finished_locus()
        return records, processor

    def test_splicegirl_read_support_is_the_junction_ratio(self):
        """ SpliceGirl reuses AD for splice-supporting reads and DP for *reference* reads. The
            ^SpliceGirl VCFSourceSettings row still reaches the rewritten VCF through its '##source', so
            VAF = AD / (AD + DP) = INFO/ALTDEDUP / (ALTDEDUP + REFDEDUP), and nothing reads DP as depth """
        records, processor = self._process()
        vcf = processor.vcf
        self.assertEqual("SpliceGirl 1.0.0.614", vcf.source)
        self.assertEqual(("AD", "DP"), (vcf.alt_depth_field, vcf.ref_depth_field))

        allele_depth_index = COHORT_GENOTYPE_HEADER.index("samples_allele_depth") \
            - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS
        read_depth_index = COHORT_GENOTYPE_HEADER.index("samples_read_depth") \
            - BulkGenotypeVCFProcessor.COHORT_GT_NUM_ADDED_FIELDS
        vaf_index = processor.cohort_gt_vaf_index
        by_splice = {percent_decode_info_value(record.INFO.get(SPLICE_INFO)): cohort_genotype
                     for record, cohort_genotype in zip(records, processor.cohort_genotypes)}

        ar_v7 = by_splice["AR-V7 splice"]
        self.assertAlmostEqual(27 / 600, float(ar_v7[vaf_index].strip("{}")), places=3)
        self.assertEqual("{27}", ar_v7[allele_depth_index])
        self.assertIsNone(ar_v7[read_depth_index])
