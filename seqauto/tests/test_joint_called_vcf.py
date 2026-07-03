"""Tests for multi-trio JointCalledVCF support (SACGF/variantgrid_sapath#343)."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone
from rest_framework.exceptions import ValidationError

from seqauto.models import (
    Aligner,
    BamFile,
    Fastq,
    JointCalledVCF,
    SampleFromSequencingSample,
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    UnalignedReads,
    SingleSampleVCF,
    VariantCaller,
    VCFFromSequencingRun,
)
from seqauto.models.models_enums import DataGeneration, PairedEnd
from seqauto.serializers.sequencing_serializers import JointCalledVCFSerializer
from snpdb.models import VCF, Sample
from upload.models import BackendVCF, UploadedFile, UploadedVCF
from upload.vcf.vcf_import import link_samples_and_vcfs_to_sequencing


def _make_sequencing_run(name="RUN_001"):
    seq_model, _ = SequencerModel.objects.get_or_create(
        model="MiSeq", data_naming_convention=DataGeneration.MISEQ)
    sequencer, _ = Sequencer.objects.get_or_create(name="MS1", sequencer_model=seq_model)
    sequencing_run = SequencingRun.objects.create(name=name, sequencer=sequencer)
    return sequencing_run


def _make_sample_sheet(sequencing_run, sample_names, sheet_hash="SHA1"):
    sample_sheet = SampleSheet.objects.create(
        sequencing_run=sequencing_run,
        path=f"/data/{sequencing_run.name}/SampleSheet.csv",
        hash=sheet_hash,
    )
    SequencingRunCurrentSampleSheet.objects.create(sequencing_run=sequencing_run,
                                                   sample_sheet=sample_sheet)
    sequencing_samples = []
    for i, name in enumerate(sample_names, start=1):
        ss = SequencingSample.objects.create(
            sample_sheet=sample_sheet,
            sample_id=name,
            sample_name=name,
            sample_number=i,
            barcode="ACGT",
        )
        sequencing_samples.append(ss)
    return sample_sheet, sequencing_samples


def _make_user():
    return User.objects.create(username="testuser")


def _make_vcf(user, name, sample_names, source_path):
    """Create VCF + Samples + UploadedFile + UploadedVCF wired together (no real file IO)."""
    vcf = VCF.objects.create(name=name, date=timezone.now(), user=user,
                             genotype_samples=len(sample_names))
    samples = []
    for sname in sample_names:
        sample = Sample.objects.create(vcf=vcf, name=sname, vcf_sample_name=sname)
        samples.append(sample)
    uploaded_file = UploadedFile.objects.create(
        path=source_path, name=name, user=user, import_source="S")  # ImportSource.SEQAUTO
    uploaded_vcf = UploadedVCF.objects.create(uploaded_file=uploaded_file, vcf=vcf)
    return vcf, samples, uploaded_vcf


def _link_samples(sequencing_samples, samples):
    """Establish SampleFromSequencingSample 1:1 links."""
    for seq_sample, sample in zip(sequencing_samples, samples):
        SampleFromSequencingSample.objects.create(sample=sample, sequencing_sample=seq_sample)


class JointCalledVCFSerializerTests(TestCase):
    """Multi-trio joint-called VCF API: idempotency and conflict handling."""

    def setUp(self):
        self.sequencing_run = _make_sequencing_run("API_RUN_001")
        self.sample_sheet, _ = _make_sample_sheet(self.sequencing_run,
                                                  ["S1", "S2", "S3"],
                                                  sheet_hash="HSH001")
        self.caller = VariantCaller.objects.create(name="gatk", version="4.0")
        self.lookup = {
            "sequencing_run": self.sequencing_run.name,
            "hash": self.sample_sheet.hash,
        }

    def _payload(self, path, caller_data=None):
        return {
            "path": path,
            "sample_sheet": self.lookup,
            "variant_caller": caller_data or {"name": self.caller.name,
                                              "version": self.caller.version,
                                              "run_params": None},
        }

    def test_idempotent_repost_same_path(self):
        s1 = JointCalledVCFSerializer(data=self._payload("/d/trio_a.vcf.gz"))
        s1.is_valid(raise_exception=True)
        s1.save()
        s2 = JointCalledVCFSerializer(data=self._payload("/d/trio_a.vcf.gz"))
        s2.is_valid(raise_exception=True)
        s2.save()
        self.assertEqual(JointCalledVCF.objects.filter(path="/d/trio_a.vcf.gz").count(), 1)
        self.assertEqual(JointCalledVCF.objects.count(), 1)

    def test_two_distinct_paths_same_sheet_and_caller(self):
        for path in ("/d/trio_a.vcf.gz", "/d/trio_b.vcf.gz"):
            s = JointCalledVCFSerializer(data=self._payload(path))
            s.is_valid(raise_exception=True)
            s.save()
        rows = JointCalledVCF.objects.filter(
            sample_sheet=self.sample_sheet, variant_caller=self.caller).order_by("path")
        self.assertEqual([r.path for r in rows], ["/d/trio_a.vcf.gz", "/d/trio_b.vcf.gz"])

    def test_conflicting_variant_caller_raises(self):
        s1 = JointCalledVCFSerializer(data=self._payload("/d/trio_a.vcf.gz"))
        s1.is_valid(raise_exception=True)
        s1.save()

        # Second post: same path, different caller version → validator rejects.
        other_caller_data = {"name": self.caller.name, "version": "5.0", "run_params": None}
        s2 = JointCalledVCFSerializer(
            data=self._payload("/d/trio_a.vcf.gz", caller_data=other_caller_data))
        s2.is_valid(raise_exception=True)
        with self.assertRaises(ValidationError):
            s2.save()
        # Original row unchanged, still referencing the original caller.
        existing = JointCalledVCF.objects.get(path="/d/trio_a.vcf.gz")
        self.assertEqual(existing.variant_caller_id, self.caller.pk)


class LinkSamplesJointCallPreservedTests(TestCase):
    """A subsequent single-sample VCF upload must not delete the joint-called VCF's grid row."""

    def setUp(self):
        self.user = _make_user()
        self.sequencing_run = _make_sequencing_run("LINK_RUN_001")
        self.sample_sheet, self.sequencing_samples = _make_sample_sheet(
            self.sequencing_run, ["S1", "S2", "S3"], sheet_hash="HSH002")
        self.caller = VariantCaller.objects.create(name="gatk", version="4.0")

    def _make_joint_called_backend_vcf(self, path, sample_names):
        vcf, samples, uploaded_vcf = _make_vcf(self.user, "trio_a", sample_names, path)
        # Link VCF samples to the matching SequencingSamples
        seq_samples_by_name = {s.sample_name: s for s in self.sequencing_samples}
        _link_samples([seq_samples_by_name[n] for n in sample_names], samples)
        joint_called_vcf = JointCalledVCF.objects.create(
            sequencing_run=self.sequencing_run,
            sample_sheet=self.sample_sheet,
            variant_caller=self.caller,
            path=path,
        )
        backend = BackendVCF.objects.create(uploaded_vcf=uploaded_vcf,
                                            joint_called_vcf=joint_called_vcf,
                                            single_sample_vcf=None)
        return backend, samples

    def _make_single_sample_backend_vcf(self, path, sample_name):
        vcf, samples, uploaded_vcf = _make_vcf(self.user, sample_name, [sample_name], path)
        # Need a BamFile to satisfy SingleSampleVCF FK. Build the chain minimally.
        seq_sample = next(s for s in self.sequencing_samples if s.sample_name == sample_name)
        fastq = Fastq.objects.create(path=f"/d/{sample_name}.fastq.gz",
                                     sequencing_sample=seq_sample,
                                     name=f"{sample_name}.fastq.gz",
                                     read=PairedEnd.R1,
                                     sequencing_run=self.sequencing_run)
        unaligned = UnalignedReads.objects.create(sequencing_sample=seq_sample,
                                                  fastq_r1=fastq)
        aligner, _ = Aligner.objects.get_or_create(name="bwa", version="0.7")
        bam_file = BamFile.objects.create(path=f"/d/{sample_name}.bam",
                                          name=f"{sample_name}.bam",
                                          sequencing_run=self.sequencing_run,
                                          unaligned_reads=unaligned,
                                          aligner=aligner)
        single_sample_vcf = SingleSampleVCF.objects.create(path=path,
                                                           sequencing_run=self.sequencing_run,
                                                           bam_file=bam_file,
                                                           variant_caller=self.caller)
        backend = BackendVCF.objects.create(uploaded_vcf=uploaded_vcf,
                                            joint_called_vcf=None,
                                            single_sample_vcf=single_sample_vcf)
        return backend, samples

    def test_joint_called_survives_single_sample_reupload(self):
        # Step 1: joint call covers S1/S2/S3.
        joint_backend, _ = self._make_joint_called_backend_vcf(
            "/d/trio_a.vcf.gz", ["S1", "S2", "S3"])
        link_samples_and_vcfs_to_sequencing(joint_backend)
        joint_row = VCFFromSequencingRun.objects.get(vcf=joint_backend.vcf)
        self.assertEqual(joint_row.variant_caller_id, self.caller.pk)

        # Step 2: a single-sample VCF for S1 arrives later.
        single_backend, _ = self._make_single_sample_backend_vcf(
            "/d/S1_single.vcf.gz", "S1")
        link_samples_and_vcfs_to_sequencing(single_backend)

        # Joint-called row must still be there — single-sample re-upload should not sweep it.
        self.assertTrue(
            VCFFromSequencingRun.objects.filter(vcf=joint_backend.vcf).exists(),
            "Joint-called VCFFromSequencingRun was deleted by the single-sample branch",
        )
        # And the single-sample row exists too.
        self.assertTrue(
            VCFFromSequencingRun.objects.filter(vcf=single_backend.vcf).exists())
        self.assertEqual(VCFFromSequencingRun.objects.filter(
            sequencing_run=self.sequencing_run).count(), 2)
