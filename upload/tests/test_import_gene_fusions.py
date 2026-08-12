"""End-to-end import of AllFusions.csv - the plan's "Done when" criteria."""
import os

import cyvcf2

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.genomics.vcf_utils import vcf_get_ref_alt_svlen_and_modification
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from genes.models import GeneFusion, FusionGeneId
from genes.tests.test_gene_fusions import GeneFusionTestCase
from snpdb.models import CohortGenotype, GenomeBuild, ImportSource, ImportStatus, Sample
from snpdb.models.models_enums import CohortGenotypeCollectionType
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
    GeneFusionCreateVCFModelTask,
    GeneFusionCreateVCFTask,
    GeneFusionInsertTask,
    _resolve_rows,
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


class TestImportGeneFusions(GeneFusionTestCase):

    def setUp(self):
        super().setUp()
        get_fake_annotation_version(GenomeBuild.get_name_or_alias("GRCh37"))
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

        # Step 1 - write the VCF the insert pipeline consumes
        create_vcf_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                    name="Create Gene Fusion Variant VCF", sort_order=0,
                                                    input_filename=ALL_FUSIONS_CSV,
                                                    output_filename=self.vcf_filename)
        self.rows_processed = GeneFusionCreateVCFTask.process_items(create_vcf_step)

        # Step 2 - the VCF/Sample/Cohort the rows hang off
        UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                  name=UploadStep.PREPROCESS_VCF_NAME, sort_order=1)
        model_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                               name="Create Data from VCF Header", sort_order=2)
        GeneFusionCreateVCFModelTask.process_items(model_step)

        # The pipeline inserts the variants between steps 2 and 3 - stand in for it, since the insert
        # path itself is the standard one and tested with the rest of the VCF pipeline
        for resolved_fusion, _row in _resolve_rows(read_all_fusions(ALL_FUSIONS_CSV)[1]):
            resolved_fusion.get_or_create()

        # Step 3 - GeneFusions and CohortGenotypes against the now-existing variants
        insert_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                name="Gene Fusion Insert", sort_order=3)
        GeneFusionInsertTask.process_items(insert_step)
        self.sample = Sample.objects.get(vcf__uploadedvcf__file_upload=self.file_upload)

    def _cohort_genotypes(self):
        collection = self.sample.vcf.cohort.cohortgenotypecollection_set.get(
            collection_type=CohortGenotypeCollectionType.UNCOMMON)
        return CohortGenotype.objects.filter(collection=collection)

    def test_all_rows_import(self):
        self.assertEqual(EXPECTED_ROWS, self.rows_processed)

    def test_creates_its_own_vcf_and_sample(self):
        vcf = self.sample.vcf
        self.assertEqual(ImportStatus.SUCCESS, vcf.import_status)
        self.assertEqual(ImportStatus.SUCCESS, self.sample.import_status)
        self.assertEqual("FusionProcessor 1.0.0.614", vcf.source)
        self.assertEqual(1, vcf.sample_set.count())

    def test_sept14_resolves_to_septin14(self):
        """ The file says SEPT14; the fusion is EGFR-SEPTIN14 """
        egfr = FusionGeneId.objects.get(symbol_str="EGFR")
        gene_fusion = GeneFusion.objects.get(anchor=egfr)
        self.assertEqual("EGFR-SEPTIN14", gene_fusion.canonical_str)

    def test_repeated_gene_pair_collapses_to_one_fusion(self):
        """ ENTPD3-RPL14 appears three times from one caller with three different 5' breakpoints """
        entpd3 = FusionGeneId.objects.get(symbol_str="ENTPD3")
        gene_fusion = GeneFusion.objects.get(anchor=entpd3, is_ordered=True)
        cohort_genotype = self._cohort_genotypes().get(variant=gene_fusion.variant)
        self.assertEqual(3, len(cohort_genotype.info["observations"]))

    def test_every_breakpoint_is_preserved(self):
        entpd3 = FusionGeneId.objects.get(symbol_str="ENTPD3")
        gene_fusion = GeneFusion.objects.get(anchor=entpd3, is_ordered=True)
        cohort_genotype = self._cohort_genotypes().get(variant=gene_fusion.variant)
        breakpoints = {o["Gene A Breakpoint"] for o in cohort_genotype.info["observations"]}
        self.assertEqual({"chr3:40442308", "chr3:40446552", "chr3:40447887"}, breakpoints)

    def test_reciprocal_pair_stays_separate(self):
        """ The file has ENTPD3->RPL14 and RPL14->ENTPD3 as separate calls """
        entpd3 = FusionGeneId.objects.get(symbol_str="ENTPD3")
        rpl14 = FusionGeneId.objects.get(symbol_str="RPL14")
        self.assertTrue(GeneFusion.objects.filter(anchor=entpd3, partner=rpl14, is_ordered=True).exists())
        self.assertTrue(GeneFusion.objects.filter(anchor=rpl14, partner=entpd3, is_ordered=True).exists())

    def test_merge_is_recorded(self):
        miv_qs = ModifiedImportedVariant.objects.filter(
            operation=ModifiedImportedVariantOperation.SHARED_LOCUS)
        self.assertTrue(miv_qs.exists(), "the 3-row ENTPD3-RPL14 merge is recorded")
        self.assertIn("3 calls merged", miv_qs.first().operation_detail)

    def test_multi_gene_partner_resolves(self):
        """ CD74 | ROS1;GOPC -> CD74-ROS1 """
        cd74 = FusionGeneId.objects.get(symbol_str="CD74")
        gene_fusion = GeneFusion.objects.get(anchor=cd74)
        self.assertEqual("CD74-ROS1", gene_fusion.canonical_str)
        cohort_genotype = self._cohort_genotypes().get(variant=gene_fusion.variant)
        self.assertEqual("ROS1;GOPC", cohort_genotype.info["observations"][0]["Gene B"],
                         "the cell is kept exactly as the caller wrote it")

    def test_clone_based_partner_gets_a_local_identity(self):
        """ PPARG/AC016683.6 -> PPARG, and RP11-458D21.5;NOTCH2NL -> NOTCH2NL """
        pparg = FusionGeneId.objects.get(symbol_str="PPARG")
        self.assertTrue(GeneFusion.objects.filter(anchor=pparg).exists())
        notch2nl = FusionGeneId.objects.get(symbol_str="NOTCH2NL")
        self.assertTrue(GeneFusion.objects.filter(anchor=notch2nl).exists())

    def test_variants_are_on_the_gene_level_contig(self):
        for cohort_genotype in self._cohort_genotypes():
            self.assertTrue(cohort_genotype.variant.is_gene_level)

    def test_one_cohort_genotype_per_fusion(self):
        cohort_genotypes = self._cohort_genotypes()
        variant_ids = {cg.variant_id for cg in cohort_genotypes}
        self.assertEqual(len(variant_ids), cohort_genotypes.count())
        self.assertLess(cohort_genotypes.count(), EXPECTED_ROWS, "some rows share a gene pair")

    def test_written_vcf_is_readable_by_the_pipeline(self):
        """ The insert pipeline reads the VCF with cyvcf2 and derives svlen the same way it does for
            any symbolic alt - END=POS is what gives gene-level variants svlen 0 """
        reader = cyvcf2.VCF(self.vcf_filename)
        records = list(reader)
        self.assertEqual(self._cohort_genotypes().count(), len(records),
                         "one record per distinct gene pair")
        for record in records:
            ref, alt, svlen, _modification = vcf_get_ref_alt_svlen_and_modification(
                record, old_variant_info=ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            self.assertEqual("N", ref)
            self.assertIsNotNone(GeneLevelSymbolicAlt.parse(alt), alt)
            self.assertEqual(0, svlen)

    def test_written_vcf_is_sorted(self):
        """ Preprocess skips vcf_clean_and_filter, which is what would otherwise catch an unsorted VCF """
        positions = [record.POS for record in cyvcf2.VCF(self.vcf_filename)]
        self.assertEqual(sorted(positions), positions)
