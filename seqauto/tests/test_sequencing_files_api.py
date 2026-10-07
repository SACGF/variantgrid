"""Tests for sequencing data sent up without FastQs (SACGF/variantgrid_sapath#357)."""
from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse
from rest_framework.exceptions import ValidationError

from seqauto.grids.sequencing_data_grids import AlignmentFileColumns
from seqauto.models import (
    QC,
    AlignmentFile,
    EnrichmentKit,
    Fastq,
    QCGeneCoverage,
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    SingleSampleVCF,
    UnalignedReads,
)
from seqauto.models.models_enums import DataGeneration
from seqauto.models.models_seqauto import SampleFromSequencingSample
from seqauto.serializers.seqauto_qc_serializers import QCSerializer
from seqauto.serializers.sequencing_serializers import (
    AlignmentFileSerializer,
    SequencingFilesBulkCreateSerializer,
)
from snpdb.fake_data import create_fake_cohort
from snpdb.models.models_enums import SampleFileType
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_vcf import SampleFilePath

SAMPLE_NAMES = ["fake_sample_1", "fake_sample_2"]


class SequencingFilesBulkCreateTests(TestCase):
    """ The bulk create API must accept a BAM against a SequencingSample with or without FastQs """

    def setUp(self):
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer, _ = Sequencer.objects.get_or_create(name="M02027", sequencer_model=sequencer_model)
        self.enrichment_kit = EnrichmentKit.objects.create(name="idt_haem")
        self.sequencing_run = SequencingRun.objects.create(name="RUN_357", sequencer=sequencer,
                                                           enrichment_kit=self.enrichment_kit)
        self.sample_sheet = SampleSheet.objects.create(sequencing_run=self.sequencing_run,
                                                       path="/data/RUN_357/SampleSheet.csv",
                                                       hash="HASH357")
        SequencingRunCurrentSampleSheet.objects.create(sequencing_run=self.sequencing_run,
                                                       sample_sheet=self.sample_sheet)
        for i, sample_name in enumerate(SAMPLE_NAMES, start=1):
            SequencingSample.objects.create(sample_sheet=self.sample_sheet,
                                            sample_id=sample_name,
                                            sample_name=sample_name,
                                            sample_number=i,
                                            barcode="ACGT")
        self.sample_sheet_lookup = {
            "sequencing_run": self.sequencing_run.name,
            "hash": self.sample_sheet.hash,
        }
        self.user = User.objects.create(username="testuser_357")

    def _sequencing_sample_lookup(self, sample_name):
        return {"sample_sheet": self.sample_sheet_lookup, "sample_name": sample_name}

    def _record(self, sample_name, fastqs=True, bam_sequencing_sample=None, alignment_extension="bam"):
        record = {
            "sample_name": sample_name,
            "bam_file": {
                "path": f"/data/RUN_357/1_BAM/{sample_name}.hg38.{alignment_extension}",
                "aligner": {"name": "BWA", "version": "0.7.18"},
            },
            "vcf_file": {
                "path": f"/data/RUN_357/2_variants/{sample_name}.gatk.hg38.vcf.gz",
                "variant_caller": {"name": "GATK", "version": "4.1.9.0", "run_params": None},
            },
        }
        if fastqs:
            record["unaligned_reads"] = {
                "fastq_r1": {"path": f"/data/RUN_357/0_fastq/{sample_name}_R1.fastq.gz"},
                "fastq_r2": {"path": f"/data/RUN_357/0_fastq/{sample_name}_R2.fastq.gz"},
            }
        if bam_sequencing_sample:
            record["bam_file"]["sequencing_sample"] = self._sequencing_sample_lookup(bam_sequencing_sample)
        return record

    def _bulk_create(self, *records):
        serializer = SequencingFilesBulkCreateSerializer(data={"sample_sheet": self.sample_sheet_lookup,
                                                               "records": list(records)})
        serializer.is_valid(raise_exception=True)
        return serializer.save()

    def _make_qc(self, sample_name):
        """ Goes through the QC API path, which resolves the sample then finds its BAM/VCF """
        qc_data = {
            "sequencing_sample": self._sequencing_sample_lookup(sample_name),
            "alignment_files": [{"path": f"/data/RUN_357/1_BAM/{sample_name}.hg38.bam"}],
            "vcf_file": {"path": f"/data/RUN_357/2_variants/{sample_name}.gatk.hg38.vcf.gz"},
            "path": f"/data/RUN_357/4_QC/exec_stats/{sample_name}_qc_summary.txt",
        }
        return QCSerializer.get_object(qc_data)

    def _alignment_file_grid_rows(self):
        url = reverse('alignment_file_datatable')
        request = RequestFactory().get(url)
        request.resolver_match = resolve(url)
        request.user = self.user
        config = AlignmentFileColumns(request)
        return config.get_initial_queryset().values(*config.value_columns())

    def test_bam_without_fastqs(self):
        sample_name = SAMPLE_NAMES[0]
        self._bulk_create(self._record(sample_name, fastqs=False))

        sequencing_sample = SequencingSample.objects.get(sample_sheet=self.sample_sheet, sample_name=sample_name)
        self.assertFalse(UnalignedReads.objects.filter(sequencing_sample=sequencing_sample).exists())
        self.assertFalse(Fastq.objects.filter(sequencing_sample=sequencing_sample).exists())

        bam_file = sequencing_sample.alignmentfile_set.get()
        self.assertIsNotNone(bam_file, "SequencingSample finds its FastQ-less BAM")
        self.assertIsNone(bam_file.unaligned_reads)
        self.assertEqual(bam_file.sequencing_run, self.sequencing_run)

        qc = self._make_qc(sample_name)
        self.assertEqual(qc.sequencing_sample, sequencing_sample)

        # QCGeneCoverage is reached from the sample via the same path the gene coverage queries use
        QCGeneCoverage.objects.create(qc=qc, sequencing_run=self.sequencing_run,
                                      path=f"/data/RUN_357/4_QC/{sample_name}.per_gene_coverage.tsv")
        coverage_qs = QCGeneCoverage.objects.filter(qc__sequencing_sample=sequencing_sample)
        self.assertEqual(coverage_qs.count(), 1)

        grid_rows = list(self._alignment_file_grid_rows())
        self.assertEqual(len(grid_rows), 1)
        self.assertEqual(grid_rows[0]["sequencing_sample__sample_sheet__sequencing_run__name"],
                         self.sequencing_run.name)

    def test_bam_with_fastqs(self):
        sample_name = SAMPLE_NAMES[0]
        self._bulk_create(self._record(sample_name))

        sequencing_sample = SequencingSample.objects.get(sample_sheet=self.sample_sheet, sample_name=sample_name)
        bam_file = sequencing_sample.alignmentfile_set.get()
        self.assertIsNotNone(bam_file)
        self.assertEqual(bam_file.unaligned_reads.sequencing_sample, sequencing_sample)
        self.assertEqual(Fastq.objects.filter(sequencing_sample=sequencing_sample).count(), 2)

        qc = self._make_qc(sample_name)
        self.assertEqual(qc.sequencing_sample, sequencing_sample)

    def test_bam_params_without_fastqs(self):
        self._bulk_create(self._record(SAMPLE_NAMES[0], fastqs=False))
        bam_file = AlignmentFile.objects.get()
        params = bam_file.get_params()
        self.assertEqual(params["aligned_pattern"], f"{SAMPLE_NAMES[0]}_S1")
        self.assertTrue(params["bam"].endswith(f"1_BAM/{SAMPLE_NAMES[0]}.hg38.bam"))

    def test_mismatched_sequencing_sample_rejected(self):
        """ FastQs from one sample with the BAM explicitly claiming another is an error """
        record = self._record(SAMPLE_NAMES[0], bam_sequencing_sample=SAMPLE_NAMES[1])
        with self.assertRaises(ValidationError) as cm:
            self._bulk_create(record)

        message = str(cm.exception)
        for sample_name in SAMPLE_NAMES:
            self.assertIn(sample_name, message)
        self.assertFalse(AlignmentFile.objects.exists())

    def _alignment_file(self, sample_name, path, aligner_version="0.7.18"):
        serializer = AlignmentFileSerializer(data={"path": path,
                                                   "sequencing_sample": self._sequencing_sample_lookup(sample_name),
                                                   "aligner": {"name": "BWA", "version": aligner_version}})
        serializer.is_valid(raise_exception=True)
        return serializer.save()

    def test_cram_and_multiple_alignment_files_reach_sample(self):
        """ A CRAM sent as 'bam_file', and an extra BAM sent after the VCF was imported, both open in IGV """
        sample_name = SAMPLE_NAMES[0]
        self._bulk_create(self._record(sample_name, alignment_extension="cram"))
        sequencing_sample = SequencingSample.objects.get(sample_sheet=self.sample_sheet, sample_name=sample_name)
        cram = sequencing_sample.alignmentfile_set.get()
        self.assertEqual(cram.file_type, SampleFileType.CRAM)

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        sample = create_fake_cohort(self.user, grch37, name="alignment_files").vcf.sample_set.first()
        SampleFromSequencingSample.objects.create(sample=sample, sequencing_sample=sequencing_sample)
        cram.link_to_samples()
        SampleFilePath.objects.create(sample=sample, file_type=SampleFileType.BED, file_path="/data/panel.bed")

        recal_bam_path = f"/data/RUN_357/1_BAM/{sample_name}.recal.bam"
        self._alignment_file(sample_name, recal_bam_path)
        self.assertCountEqual(sample.get_alignment_files(), [cram.path, recal_bam_path])

    def test_repost_updates_alignment_file(self):
        sample_name = SAMPLE_NAMES[0]
        path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.bam"
        self._alignment_file(sample_name, path)
        alignment_file = self._alignment_file(sample_name, path, aligner_version="0.7.19")
        self.assertEqual(AlignmentFile.objects.get().pk, alignment_file.pk)
        self.assertEqual(alignment_file.aligner.version, "0.7.19")

    def test_single_qc_with_multiple_alignment_files(self):
        sample_name = SAMPLE_NAMES[0]
        self._bulk_create(self._record(sample_name))
        self._alignment_file(sample_name, f"/data/RUN_357/1_BAM/{sample_name}.recal.bam")
        qc = self._make_qc(sample_name)
        sequencing_sample = SequencingSample.objects.get(sample_sheet=self.sample_sheet, sample_name=sample_name)
        self.assertEqual(sequencing_sample.get_single_qc(), qc)

    def _alignment_files_record(self, sample_name, *alignment_files):
        record = self._record(sample_name, fastqs=False)
        record.pop("bam_file")
        record["alignment_files"] = list(alignment_files)
        return record

    def _qc_from_api(self, qc_data):
        qc_serializer = QCSerializer(data=qc_data)
        qc_serializer.is_valid(raise_exception=True)
        return QCSerializer.get_object(qc_serializer.validated_data)

    def _qc_data(self, sample_name, **kwargs):
        return {
            "sequencing_sample": self._sequencing_sample_lookup(sample_name),
            "vcf_file": {"path": f"/data/RUN_357/2_variants/{sample_name}.gatk.hg38.vcf.gz"},
            **kwargs,
        }

    def test_vcf_and_qc_link_every_alignment_file(self):
        sample_name = SAMPLE_NAMES[0]
        aligner = {"name": "BWA", "version": "0.7.18"}
        bam_path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.bam"
        cram_path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.cram"
        no_extension_path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.alignment"
        alignment_files = [{"path": bam_path, "aligner": aligner},
                           {"path": cram_path, "aligner": aligner},
                           {"path": no_extension_path, "aligner": aligner, "file_type": "cram"}]
        self._bulk_create(self._alignment_files_record(sample_name, *alignment_files))
        file_types = dict(AlignmentFile.objects.values_list("path", "file_type"))
        self.assertEqual(file_types, {bam_path: SampleFileType.BAM,
                                      cram_path: SampleFileType.CRAM,
                                      no_extension_path: SampleFileType.CRAM})

        # Order doesn't key the VCF - a re-post in another order is the same record
        self._bulk_create(self._alignment_files_record(sample_name, *reversed(alignment_files)))
        vcf_file = SingleSampleVCF.objects.get()
        all_paths = {bam_path, cram_path, no_extension_path}
        self.assertEqual(set(vcf_file.alignment_files.values_list("path", flat=True)), all_paths)

        # QC sent without alignment files takes the VCF's, a re-post naming some replaces them
        qc = self._qc_from_api(self._qc_data(sample_name))
        self.assertEqual(set(qc.alignment_files.values_list("path", flat=True)), all_paths)
        reposted_qc = self._qc_from_api(self._qc_data(sample_name, alignment_files=[{"path": cram_path}]))
        self.assertEqual(reposted_qc.pk, qc.pk)
        self.assertEqual(list(reposted_qc.alignment_files.values_list("path", flat=True)), [cram_path])

    def test_bam_file_and_alignment_files_merged(self):
        sample_name = SAMPLE_NAMES[0]
        record = self._record(sample_name)
        cram_path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.cram"
        record["alignment_files"] = [{"path": cram_path, "aligner": {"name": "BWA", "version": "0.7.18"}}]
        result = self._bulk_create(record)

        sequencing_sample = SequencingSample.objects.get(sample_sheet=self.sample_sheet, sample_name=sample_name)
        self.assertTrue(all(af.unaligned_reads for af in sequencing_sample.alignmentfile_set.all()))
        vcf_file = result["records"][0]["vcf_file"]
        self.assertEqual(set(vcf_file.alignment_files.values_list("path", flat=True)),
                         {record["bam_file"]["path"], cram_path})

    def test_alignment_files_required(self):
        with self.assertRaises(ValidationError):
            self._bulk_create(self._alignment_files_record(SAMPLE_NAMES[0]))

    def test_vcf_alignment_files_from_one_sequencing_sample(self):
        aligner = {"name": "BWA", "version": "0.7.18"}
        alignment_files = [{"path": f"/data/RUN_357/1_BAM/{sample_name}.bam", "aligner": aligner,
                            "sequencing_sample": self._sequencing_sample_lookup(sample_name)}
                           for sample_name in SAMPLE_NAMES]
        with self.assertRaises(ValidationError):
            self._bulk_create(self._alignment_files_record(SAMPLE_NAMES[0], *alignment_files))

    def test_aligner_only_needed_for_new_alignment_file(self):
        sample_name = SAMPLE_NAMES[0]
        path = f"/data/RUN_357/1_BAM/{sample_name}.hg38.bam"
        with self.assertRaises(ValidationError):
            self._bulk_create(self._alignment_files_record(sample_name, {"path": path}))

        self._bulk_create(self._record(sample_name, fastqs=False))
        self._bulk_create(self._alignment_files_record(sample_name, {"path": path}))
        self.assertEqual(AlignmentFile.objects.get().aligner.version, "0.7.18")

    def test_qc_bam_file_still_accepted(self):
        sample_name = SAMPLE_NAMES[0]
        self._bulk_create(self._record(sample_name))
        bam_file = {"path": f"/data/RUN_357/1_BAM/{sample_name}.hg38.bam"}
        qc = self._qc_from_api(self._qc_data(sample_name, bam_file=bam_file))
        self.assertEqual(QCSerializer(qc).data["bam_file"], bam_file)

        with self.assertRaises(AlignmentFile.DoesNotExist):
            self._qc_from_api(self._qc_data(sample_name, alignment_files=[{"path": "/data/other.bam"}]))
