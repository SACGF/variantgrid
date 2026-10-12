"""Sequencing run and sample sheet maintenance: re-posting a sheet, the run-created signal, the 20x gene coverage
count and the superuser-only run actions (SACGF/variantgrid_private#3913)."""
from django.contrib.auth.models import User
from django.core.cache import cache
from django.urls import reverse
from rest_framework.test import APITestCase

from genes.models import GeneCoverageCanonicalTranscript, GeneCoverageCollection, GeneSymbol
from seqauto.models import (
    QC,
    Aligner,
    AlignmentFile,
    EnrichmentKit,
    Experiment,
    QCGeneCoverage,
    SampleSheet,
    SequencingSample,
    SingleSampleVCF,
    VariantCaller,
    get_20x_gene_coverage,
)
from seqauto.signals.signals_list import sequencing_run_created_signal
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run
from snpdb.models import DataState, GenomeBuild

SAMPLE_NAMES = ["sample_a", "sample_b"]


def _sheet_payload(sequencing_run, sample_sheet, sample_names, enrichment_kit):
    kit = {"name": enrichment_kit.name, "version": enrichment_kit.version}
    return {
        "path": sample_sheet.path,
        "sequencing_run": sequencing_run.pk,
        "hash": sample_sheet.hash,
        "sequencingsample_set": [
            {"sample_id": name, "sample_name": name, "sample_number": i, "lane": 1, "barcode": "ACGT",
             "enrichment_kit": kit}
            for i, name in enumerate(sample_names, start=1)
        ],
    }


class SampleSheetUpdateTest(APITestCase):
    """ A PUT to a sample sheet updates its rows in place and keeps the data linked to them """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_superuser(username="sample_sheet_update_user")
        cls.sequencing_run = make_sequencing_run("SHEET_UPDATE_RUN")
        cls.sample_sheet, cls.sequencing_samples = make_sample_sheet(cls.sequencing_run, SAMPLE_NAMES)
        cls.enrichment_kit = EnrichmentKit.objects.create(name="sheet_update_kit")
        aligner = Aligner.objects.get_or_create(name="bwa", version="0.7")[0]
        cls.alignment_file = AlignmentFile.objects.create(path="/data/SHEET_UPDATE_RUN/1_BAM/sample_a.bam", name="sample_a.bam",
                                              sequencing_run=cls.sequencing_run,
                                              sequencing_sample=cls.sequencing_samples[0], aligner=aligner)

    def test_put_keeps_rows_and_their_data(self):
        self.client.force_authenticate(user=self.user)
        payload = _sheet_payload(self.sequencing_run, self.sample_sheet, SAMPLE_NAMES, self.enrichment_kit)
        payload["sequencingsample_set"][0]["failed"] = True
        url = reverse("api_sample_sheet-detail", kwargs={"pk": self.sample_sheet.pk})
        response = self.client.put(url, payload, format="json")
        self.assertEqual(response.status_code, 200, response.content)

        sample_a = SequencingSample.objects.get(pk=self.sequencing_samples[0].pk)
        self.assertTrue(sample_a.failed)
        self.assertTrue(AlignmentFile.objects.filter(pk=self.alignment_file.pk).exists())
        self.assertEqual(self.sample_sheet.sequencingsample_set.count(), len(SAMPLE_NAMES))


class SequencingRunCreatedSignalTest(APITestCase):
    """ The API sends sequencing_run_created_signal once, when the run is new """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_superuser(username="run_created_signal_user")
        cls.sequencer = make_sequencing_run("RUN_CREATED_SIGNAL_EXISTING").sequencer
        cls.experiment = Experiment.objects.create(name="run_created_signal_experiment")
        cls.enrichment_kit = EnrichmentKit.objects.create(name="run_created_signal_kit")

    def test_signal_sent_for_new_run_only(self):
        received = []

        def receiver(sender, sequencing_run, **kwargs):  # pylint: disable=unused-argument
            received.append(sequencing_run)

        sequencing_run_created_signal.connect(receiver)
        self.addCleanup(sequencing_run_created_signal.disconnect, receiver)

        self.client.force_authenticate(user=self.user)
        payload = {"name": "RUN_CREATED_SIGNAL_NEW", "path": "/data/RUN_CREATED_SIGNAL_NEW",
                   "sequencer": self.sequencer.pk, "experiment": self.experiment.pk,
                   "enrichment_kit": {"name": self.enrichment_kit.name, "version": self.enrichment_kit.version}}
        for _ in range(2):
            response = self.client.post(reverse("api_sequencing_run-list"), payload, format="json")
            self.assertEqual(response.status_code, 201, response.content)

        self.assertEqual([sr.pk for sr in received], ["RUN_CREATED_SIGNAL_NEW"])


class Gene20xCoverageCountTest(APITestCase):
    """ get_20x_gene_coverage counts each current-sheet collection once, across cached calls """
    GENE_SYMBOL = "RUNX1"
    MIN_COVERAGE = 100

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch38()
        cls.gene_symbol = GeneSymbol.objects.get_or_create(symbol=cls.GENE_SYMBOL)[0]
        cls.aligner = Aligner.objects.get_or_create(name="bwa", version="0.7")[0]
        cls.variant_caller = VariantCaller.objects.get_or_create(name="gatk", version="4")[0]

    def setUp(self):
        cache.delete(f"20x_gene_coverage_{self.GENE_SYMBOL}_cov_{self.MIN_COVERAGE}")

    def _collection_for(self, sequencing_run, sample_sheet, sample_name, percent_20x=100):
        sequencing_sample = SequencingSample.objects.create(sample_sheet=sample_sheet, sample_id=sample_name,
                                                            sample_name=sample_name, sample_number=1, barcode="ACGT")
        alignment_file = AlignmentFile.objects.create(path=f"/data/{sequencing_run.name}/{sample_name}.bam", name=sample_name,
                                          sequencing_run=sequencing_run, sequencing_sample=sequencing_sample,
                                          aligner=self.aligner)
        vcf_file = SingleSampleVCF.objects.create(path=f"/data/{sequencing_run.name}/{sample_name}.vcf",
                                                  sequencing_run=sequencing_run, alignment_file=alignment_file,
                                                  variant_caller=self.variant_caller)
        qc = QC.objects.create(path=f"/data/{sequencing_run.name}/{sample_name}_qc.txt", sequencing_run=sequencing_run,
                               alignment_file=alignment_file, vcf_file=vcf_file)
        collection = GeneCoverageCollection.objects.create(path=f"/data/{sequencing_run.name}/{sample_name}.cov.tsv",
                                                           data_state=DataState.COMPLETE,
                                                           genome_build=self.genome_build)
        QCGeneCoverage.objects.create(qc=qc, sequencing_run=sequencing_run,
                                      path=f"/data/{sequencing_run.name}/{sample_name}.cov.tsv",
                                      gene_coverage_collection=collection)
        GeneCoverageCanonicalTranscript.objects.create(gene_coverage_collection=collection,
                                                       gene_symbol=self.gene_symbol,
                                                       original_gene_symbol=self.GENE_SYMBOL, original_transcript="",
                                                       min=20, mean=50.0, std_dev=1.0, percent_20x=percent_20x)
        return collection

    def test_counts_current_sheet_collections_once(self):
        run_a = make_sequencing_run("COVERAGE_RUN_A")
        sheet_a, _ = make_sample_sheet(run_a, [])
        self._collection_for(run_a, sheet_a, "a1")
        self.assertEqual(get_20x_gene_coverage(self.GENE_SYMBOL, self.MIN_COVERAGE), 1)

        # A collection under a sheet that is no longer current is not counted
        run_b = make_sequencing_run("COVERAGE_RUN_B")
        old_sheet_b = SampleSheet.objects.create(sequencing_run=run_b, hash="OLD", path="/data/COVERAGE_RUN_B/old.csv")
        make_sample_sheet(run_b, [], sheet_hash="NEW")
        self._collection_for(run_b, old_sheet_b, "b_old")
        self.assertEqual(get_20x_gene_coverage(self.GENE_SYMBOL, self.MIN_COVERAGE), 1)

        # A new current collection adds exactly one, without recounting the cached maximum
        run_c = make_sequencing_run("COVERAGE_RUN_C")
        sheet_c, _ = make_sample_sheet(run_c, [])
        self._collection_for(run_c, sheet_c, "c1")
        self.assertEqual(get_20x_gene_coverage(self.GENE_SYMBOL, self.MIN_COVERAGE), 2)
        self.assertEqual(get_20x_gene_coverage(self.GENE_SYMBOL, self.MIN_COVERAGE), 2)


class SequencingRunActionPermissionTest(APITestCase):
    """ Relinking a run's data and reloading its experiment name are superuser actions """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="run_action_user")
        cls.sequencing_run = make_sequencing_run("RUN_ACTION_PERMISSION")
        make_sample_sheet(cls.sequencing_run, SAMPLE_NAMES)

    def test_non_superuser_gets_403(self):
        self.client.force_login(self.user)
        for url_name in ("assign_data_to_current_sample_sheet", "reload_experiment_name"):
            url = reverse(url_name, kwargs={"sequencing_run_id": self.sequencing_run.pk})
            self.assertEqual(self.client.post(url).status_code, 403, url_name)
