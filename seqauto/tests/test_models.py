import logging
import os

from django.conf import settings
from django.test import TestCase

from genes.canonical_transcripts.canonical_transcript_manager import CanonicalTranscriptManager
from genes.canonical_transcripts.create_canonical_transcripts import create_canonical_transcript_collection
from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneCoverageCollection, TranscriptVersion
from seqauto.models import SequencerModel, Sequencer, SequencingRun, SampleSheet, \
    SequencingRunCurrentSampleSheet, SequencingSample, Fastq, UnalignedReads, Aligner, \
    BamFile, VCFFile, QC, VariantCaller, EnrichmentKit, QCGeneCoverage
from seqauto.models.models_enums import DataGeneration
from seqauto.tasks.gold_summary_tasks import calculate_gold_summary
from snpdb.models import Manufacturer, GenomeBuild, DataState


class TestSeqAutoModels(TestCase):
    SEQAUTO_TEST_BASE_DIR = os.path.join(settings.BASE_DIR, "seqauto", "test_data")
    TEST_DATA = os.path.join(SEQAUTO_TEST_BASE_DIR, "clinical_hg38")
    SEQUENCING_RUN_CAPTURE = "Exome_20_022_200920_NB501009_0410_AHNLYFBGXG"
    CANONICAL_TRANSCRIPTS = os.path.join(SEQAUTO_TEST_BASE_DIR, "reference_data", "canonical", "fake_kit.GeneTable.tsv")

    @staticmethod
    def _create_sequencing_run(sequencing_run_name, enrichment_kit=None):
        sequencer_model, _ = SequencerModel.objects.get_or_create(model='NextSeq',
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer, _ = Sequencer.objects.get_or_create(name="my sequencer",
                                                       sequencer_model=sequencer_model)
        sequencing_run, _ = SequencingRun.objects.get_or_create(name=sequencing_run_name,
                                                                sequencer=sequencer,
                                                                enrichment_kit=enrichment_kit)
        return sequencing_run

    @staticmethod
    def _create_sequencing_qc(sequencing_run, sample_name, enrichment_kit=None):
        """ Create all the way down to QC """

        sample_sheet, _ = SampleSheet.objects.get_or_create(sequencing_run=sequencing_run)
        SequencingRunCurrentSampleSheet.objects.get_or_create(sequencing_run=sequencing_run,
                                                              sample_sheet=sample_sheet)
        sequencing_sample, _ = SequencingSample.objects.get_or_create(sample_sheet=sample_sheet,
                                                                      sample_name=sample_name,
                                                                      sample_number=1,
                                                                      enrichment_kit=enrichment_kit)
        fastq, _ = Fastq.objects.get_or_create(sequencing_sample=sequencing_sample)
        unaligned_reads, _ = UnalignedReads.objects.get_or_create(sequencing_sample=sequencing_sample,
                                                                  fastq_r1=fastq)
        aligner, _ = Aligner.objects.get_or_create(name="Fake Aligner")
        bam_file, _ = BamFile.objects.get_or_create(unaligned_reads=unaligned_reads,
                                                    aligner=aligner)
        variant_caller, _ = VariantCaller.objects.get_or_create(name="Fake Caller")
        vcf_file, _ = VCFFile.objects.get_or_create(bam_file=bam_file,
                                                    variant_caller=variant_caller)
        path = os.path.join(TestSeqAutoModels.TEST_DATA,
                            "idt_exome/Exome_20_022_200920_NB501009_0410_AHNLYFBGXG/4_QC/exec_stats/hiseq_sample1_stats.txt")
        qc, _ = QC.objects.get_or_create(bam_file=bam_file, vcf_file=vcf_file,
                                         defaults={"path": path})
        return qc

    def _create_enrichment_kit(self, name="fake_kit"):
        manufacturer, _ = Manufacturer.objects.get_or_create(name="Agilluminent")
        enrichment_kit, _ = EnrichmentKit.objects.get_or_create(name=name,
                                                                version=1,
                                                                manufacturer=manufacturer)

        if enrichment_kit.canonical_transcript_collection is None:
            enrichment_kit.canonical_transcript_collection = create_canonical_transcript_collection(self.genome_build,
                                                                                                    self.genome_build.annotation_consortium,
                                                                                                    TestSeqAutoModels.CANONICAL_TRANSCRIPTS,
                                                                                                    gene_matcher=GeneSymbolMatcher())
        return enrichment_kit

    @staticmethod
    def _create_qc_gene_coverage(qc, gene_filename, genome_build: GenomeBuild):
        gcc, _ = GeneCoverageCollection.objects.get_or_create(path=gene_filename, genome_build=genome_build)
        qcgc_defaults = {"path": gene_filename,
                         "data_state": DataState.COMPLETE,
                         "gene_coverage_collection": gcc}
        qcgc, _ = QCGeneCoverage.objects.get_or_create(qc=qc, defaults=qcgc_defaults)
        return qcgc

    def setUp(self):
        EnrichmentKit.objects.all().delete()
        SequencingRun.objects.all().delete()
        self.genome_build = GenomeBuild.grch37()
        self.gene_matcher = GeneSymbolMatcher()
        self.canonical_transcript_manager = CanonicalTranscriptManager(use_system_default=False)
        self.transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id(self.genome_build, self.genome_build.annotation_consortium)

    def test_load_qc_exec_summary(self):
        enrichment_kit = self._create_enrichment_kit("idt_exome")
        sequencing_run = TestSeqAutoModels._create_sequencing_run(self.SEQUENCING_RUN_CAPTURE,
                                                                  enrichment_kit=enrichment_kit)
        qc = TestSeqAutoModels._create_sequencing_qc(sequencing_run,
                                                     "hiseq_sample1",
                                                     enrichment_kit=enrichment_kit)
        qc.load_from_file(None)
        qc_exec_summary = qc.qcexecsummary_set.order_by("pk").last()
        self.assertAlmostEqual(qc_exec_summary.mean_coverage_across_genes, 262.84, places=2)

    def test_gene_coverage_capture(self):
        """ Capture uses headers: ["% bases >20x", "% bases <10x"] """

        CAPTURE_GENES = os.path.join(self.TEST_DATA, f"idt_exome/{self.SEQUENCING_RUN_CAPTURE}/4_QC/bam_stats/samples/hiseq_sample1.per_gene_coverage.tsv.gz")
        enrichment_kit = self._create_enrichment_kit()
        sequencing_run = TestSeqAutoModels._create_sequencing_run(self.SEQUENCING_RUN_CAPTURE,
                                                                  enrichment_kit=enrichment_kit)
        qc = TestSeqAutoModels._create_sequencing_qc(sequencing_run,
                                                     "hiseq_sample1",
                                                     enrichment_kit=enrichment_kit)
        qcgc = self._create_qc_gene_coverage(qc, CAPTURE_GENES, self.genome_build)
        qcgc.load_from_file(None,
                            gene_matcher=self.gene_matcher,
                            canonical_transcript_manager=self.canonical_transcript_manager,
                            transcript_versions_by_id=self.transcript_versions_by_id)
        gcc = qcgc.gene_coverage_collection
        num_genes_covered = gcc.genecoverage_set.count()
        num_canonical_genes = gcc.genecoveragecanonicaltranscript_set.count()
        logging.info("genes: %d, canonical: %d", num_genes_covered, num_canonical_genes)

    def test_deleted_gene_coverage(self):
        """ #1619 - Ensure missing file sets data_state=DELETED and doesn't throw error. """
        MISSING_GENES_FILE = os.path.join(self.TEST_DATA, "no_such_file.txt")
        enrichment_kit = self._create_enrichment_kit()
        sequencing_run = TestSeqAutoModels._create_sequencing_run(self.SEQUENCING_RUN_CAPTURE,
                                                                  enrichment_kit=enrichment_kit)
        qc = TestSeqAutoModels._create_sequencing_qc(sequencing_run,
                                                     "blah_blah",
                                                     enrichment_kit=enrichment_kit)
        qcgc = self._create_qc_gene_coverage(qc, MISSING_GENES_FILE, self.genome_build)
        qcgc.load_from_file(None,
                            gene_matcher=self.gene_matcher,
                            canonical_transcript_manager=self.canonical_transcript_manager,
                            transcript_versions_by_id=self.transcript_versions_by_id)
        msg = f"Missing file should set data_state to DELETED was {qcgc.get_data_state_display()}"
        self.assertEqual(qcgc.data_state, DataState.DELETED, msg)

    def test_gold_coverage(self):
        """ This test doesn't work very well as there's no genes/Transcripts/UCSC aliases etc """

        SEQUENCING_RUN_NAME = "Exome_20_022_200920_NB501009_0410_AHNLYFBGXG"
        GENE_FILES_PATTERN = os.path.join(self.TEST_DATA, "idt_exome", SEQUENCING_RUN_NAME, "4_QC", "bam_stats",
                                          "samples", "%s.per_gene_coverage.tsv.gz")

        enrichment_kit = self._create_enrichment_kit()
        sequencing_run = TestSeqAutoModels._create_sequencing_run(SEQUENCING_RUN_NAME,
                                                                  enrichment_kit=enrichment_kit)
        sequencing_run.gold_standard = True
        sequencing_run.save()

        SAMPLE_NAMES = ["hiseq_sample1", "hiseq_sample2"]
        for sample_name in SAMPLE_NAMES:
            qc = TestSeqAutoModels._create_sequencing_qc(sequencing_run,
                                                         sample_name,
                                                         enrichment_kit=enrichment_kit)

            gene_filename = GENE_FILES_PATTERN % sample_name
            logging.info(gene_filename)
            qcgc = self._create_qc_gene_coverage(qc, gene_filename, self.genome_build)
            print(f"created coverage for qc: {qc} - data state: {qcgc.get_data_state_display()} path: {qc.path}")
            qcgc.load_from_file(None,
                                gene_matcher=self.gene_matcher,
                                canonical_transcript_manager=self.canonical_transcript_manager,
                                transcript_versions_by_id=self.transcript_versions_by_id)

        calculate_gold_summary(enrichment_kit.pk)

        num_gold_samples = enrichment_kit.goldreference.goldgenecoveragecollection_set.count()
        sample_names = ', '.join(SAMPLE_NAMES)
        msg = f"Created gold samples for each of {sample_names}"
        self.assertEqual(num_gold_samples, len(SAMPLE_NAMES), msg)
