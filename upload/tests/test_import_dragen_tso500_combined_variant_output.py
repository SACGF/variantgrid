"""Import of CombinedVariantOutput.tsv - the splice VCF the loader writes from it."""
import os

import cyvcf2
import simplejson
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from genes.gene_splice import coordinate_label
from genes.models import HGNC, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from library.genomics.vcf_writer import percent_decode_info_value
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from snpdb.models import VCF, GenomeBuild, ImportSource, VCFSourceSettings
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline, UploadStep
from upload.tasks.import_dragen_tso500_combined_variant_output_task import (
    ALT_READS_FORMAT,
    REF_READS_FORMAT,
    DragenTSO500CombinedVariantOutputCreateVCFTask,
)
from upload.tso500.dragen_combined_variant_output_parser import (
    AFFECTED_EXON,
    BREAKPOINT_1,
    RNA_SAMPLE_ID,
    SPLICE_INFO,
    SPLICE_OBSERVATION_INFO,
    can_process_file,
    get_analysis_details,
    get_splice_rows,
    read_combined_variant_output,
)

TSO500_PAIR_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                               "ExampleSample_2600000001")
COMBINED_VARIANT_OUTPUT = os.path.join(TSO500_PAIR_DIR,
                                       "ExampleSample_2600000001_CombinedVariantOutput.tsv")
EXPECTED_SPLICE_ROWS = 3


class TestCombinedVariantOutputParser(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.sections = read_combined_variant_output(COMBINED_VARIANT_OUTPUT)

    def test_analysis_details_are_key_values(self):
        details = get_analysis_details(self.sections)
        self.assertEqual("ExampleSample_RNA_2600000001B", details[RNA_SAMPLE_ID])
        self.assertEqual("2.1.1", details["Module Version"])

    def test_reads_every_splice_row(self):
        rows = get_splice_rows(self.sections)
        self.assertEqual(EXPECTED_SPLICE_ROWS, len(rows))
        self.assertEqual("chrX:66905968", rows[0][BREAKPOINT_1])

    def test_blank_affected_exon_is_none(self):
        """ AR-V7's 3' side is a cryptic exon, so the caller leaves the column empty """
        rows = {r["Gene"]: r for r in get_splice_rows(self.sections)}
        self.assertIsNone(rows["AR"][AFFECTED_EXON])
        self.assertEqual("2-7", rows["EGFR"][AFFECTED_EXON])

    def test_na_section_is_empty(self):
        """ '[Sequencing Run Details]' is a lone NA - a section with nothing in it """
        self.assertEqual([], self.sections["Sequencing Run Details"].lines)

    def test_section_it_has_no_name_for_is_still_read(self):
        """ 2.6 renames '[Exon-Level CNVs]' and adds sections of its own, so an unknown one is data """
        rows = self.sections["Exon-Level CNVs"].rows
        self.assertEqual(["BRCA1", "BRCA2"], [r["Gene"] for r in rows])

    def test_tab_padding_is_dropped(self):
        """ Every line is padded to the widest section, 11 columns """
        self.assertEqual(["Total TMB", "7.1"], self.sections["TMB"].lines[0])

    def test_claims_the_file_over_a_gene_list(self):
        factories = {type(f).__name__: f for f in get_import_task_factories()}
        combined = factories["DragenTSO500CombinedVariantOutputImportTaskFactory"]
        gene_list = factories["GeneListImportTaskFactory"]
        user = User.objects.get_or_create(username='testuser')[0]
        self.assertGreater(combined.get_processing_ability(user, COMBINED_VARIANT_OUTPUT, "tsv"),
                           gene_list.get_processing_ability(user, COMBINED_VARIANT_OUTPUT, "tsv"))

    def test_does_not_claim_other_tsvs(self):
        """ A gene table is a tsv of gene symbols - only the banner line says this is a CVO """
        self.assertFalse(can_process_file(os.path.join(settings.BASE_DIR, "seqauto", "test_data",
                                                       "reference_data", "canonical",
                                                       "fake_kit.GeneTable.tsv")))


class TestSpliceVariantVCF(TestCase):
    """ The VCF the loader writes - what the standard insert pipeline then consumes """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(644, "AR"), (3236, "EGFR"), (7029, "MET")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk

    def setUp(self):
        self.user = User.objects.get_or_create(username='testuser')[0]
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "splice_variants.vcf")
        self.records = self._process()

    def _process(self) -> list:
        file_upload = FileUpload.objects.create(path=COMBINED_VARIANT_OUTPUT,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=self.user,
                                                name="ExampleSample_2600000001_CombinedVariantOutput.tsv",
                                                file_type=UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT,
                                                metadata={"genome_build": "GRCh37"})
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                name="Create Splice Variant VCF", sort_order=0,
                                                input_filename=COMBINED_VARIANT_OUTPUT,
                                                output_filename=self.vcf_filename)
        rows = DragenTSO500CombinedVariantOutputCreateVCFTask.process_items(upload_step)
        self.assertEqual(EXPECTED_SPLICE_ROWS, rows)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        return list(self.reader)

    def test_header_names_the_rna_sample_and_the_module(self):
        """ The splice caller runs on the RNA arm, and '##source' is what VCFSourceSettings match """
        self.assertEqual(["ExampleSample_RNA_2600000001B"], self.reader.samples)
        self.assertIn("DRAGEN TSO500 CombinedVariantOutput 2.1.1", self.reader.raw_header)

    def test_records_are_gene_level_splice_alts(self):
        self.assertEqual(EXPECTED_SPLICE_ROWS, len(self.records))
        for record in self.records:
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            kind, namespace, _gene_id, label = GeneLevelSymbolicAlt.parse(record.ALT[0])
            self.assertEqual(GeneLevelSymbolicAlt.SPLICE, kind)
            self.assertEqual(GeneIdNamespace.HGNC, namespace)
            self.assertTrue(label)

    def test_seeded_junctions_get_their_label(self):
        alts = {record.ALT[0] for record in self.records}
        self.assertEqual({f"<SPLICE:HGNC:{self.hgnc_ids['AR']}:V7>",
                          f"<SPLICE:HGNC:{self.hgnc_ids['EGFR']}:vIII>",
                          f"<SPLICE:HGNC:{self.hgnc_ids['MET']}:ex14skip>"}, alts)

    def _by_splice(self) -> dict:
        """ INFO values are stored as the VCF wrote them - a space is percent encoded, so decode """
        return {percent_decode_info_value(record.INFO.get(SPLICE_INFO)): record
                for record in self.records}

    def test_position_is_the_gene(self):
        self.assertEqual(self.hgnc_ids["AR"], self._by_splice()["AR V7"].POS)

    def test_read_support_is_the_junction_and_the_reference_transcript(self):
        """ A caller asserts the junction is present, so there is no GT - the sample column holds
            how many reads crossed it and how many crossed the reference transcript """
        support = {}
        for splice, record in self._by_splice().items():
            self.assertEqual([ALT_READS_FORMAT, REF_READS_FORMAT], record.FORMAT)
            support[splice] = (int(record.format(ALT_READS_FORMAT).flatten()[0]),
                               int(record.format(REF_READS_FORMAT).flatten()[0]))
        self.assertEqual({"AR V7": (27, 573), "EGFR vIII": (64, 1), "MET ex14skip": (91, 1)}, support)

    def test_the_callers_row_rides_along_in_info(self):
        encoded = self._by_splice()["MET ex14skip"].INFO.get(SPLICE_OBSERVATION_INFO)
        observation = simplejson.loads(percent_decode_info_value(encoded))
        self.assertEqual("chr7:116411708", observation[BREAKPOINT_1])
        self.assertEqual("14", observation[AFFECTED_EXON])

    def test_unnamed_junction_is_labelled_with_its_coordinates(self):
        """ A junction no SpliceEvent names still imports - the label reads as raw coordinates,
            which is the prompt to add a row """
        splice_event = SpliceEvent.objects.get(genome_build=self.genome_build, label="V7")
        contig = splice_event.contig
        splice_event.delete()

        records = self._process()
        expected = coordinate_label(contig, 66905968, 66914514)
        self.assertIn(f"<SPLICE:HGNC:{self.hgnc_ids['AR']}:{expected}>",
                      {record.ALT[0] for record in records})

    def test_source_settings_bind_read_support_as_depth(self):
        """ The ^DRAGEN TSO500 CombinedVariantOutput row makes ALT_READS/REF_READS the depths the
            sample node and VAF use, and with no GT in the header nothing binds a genotype """
        vcf = VCF(source="DRAGEN TSO500 CombinedVariantOutput 2.1.1", genotype_samples=1)
        for vss in VCFSourceSettings.get_for_source(vcf.source):
            vss.apply_sample_field_overrides(vcf)
        self.assertEqual(ALT_READS_FORMAT, vcf.alt_depth_field)
        self.assertEqual(REF_READS_FORMAT, vcf.ref_depth_field)
        self.assertTrue(vcf.has_depth)
        self.assertFalse(vcf.has_genotype)
