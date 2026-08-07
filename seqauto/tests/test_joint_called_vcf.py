"""Tests for multi-trio JointCalledVCF support (SACGF/variantgrid_sapath#343) and
   joint calls whose samples span sequencing runs (SACGF/variantgrid_sapath#415)."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone
from rest_framework.exceptions import ValidationError

from seqauto.models import (
    Aligner,
    BamFile,
    EnrichmentKit,
    Fastq,
    JointCalledVCF,
    SampleFromSequencingSample,
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    SingleSampleVCF,
    UnalignedReads,
    VariantCaller,
    VCFFromSequencingRun,
)
from seqauto.models.models_enums import DataGeneration, PairedEnd
from seqauto.sequencing_files.sample_sheet import current_sample_sheet_changed
from seqauto.serializers.sequencing_serializers import JointCalledVCFSerializer
from snpdb.models import VCF, Sample
from upload.models import (
    BackendVCF,
    FileUpload,
    SimpleVCFImportInfo,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
)
from upload.vcf.vcf_import import MIXED_ENRICHMENT_KIT_MESSAGE, link_samples_and_vcfs_to_sequencing


def _make_sequencing_run(name="RUN_001", enrichment_kit=None):
    seq_model, _ = SequencerModel.objects.get_or_create(
        model="MiSeq", data_naming_convention=DataGeneration.MISEQ)
    sequencer, _ = Sequencer.objects.get_or_create(name="MS1", sequencer_model=seq_model)
    sequencing_run = SequencingRun.objects.create(name=name, sequencer=sequencer,
                                                  enrichment_kit=enrichment_kit)
    return sequencing_run


def _make_sample_sheet(sequencing_run, sample_names, sheet_hash="SHA1", set_current=True):
    sample_sheet = SampleSheet.objects.create(
        sequencing_run=sequencing_run,
        path=f"/data/{sequencing_run.name}/SampleSheet_{sheet_hash}.csv",
        hash=sheet_hash,
    )
    if set_current:
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
    """Create VCF + Samples + FileUpload + UploadedVCF wired together (no real file IO)."""
    vcf = VCF.objects.create(name=name, date=timezone.now(), user=user,
                             genotype_samples=len(sample_names))
    samples = []
    for sname in sample_names:
        sample = Sample.objects.create(vcf=vcf, name=sname, vcf_sample_name=sname)
        samples.append(sample)
    file_upload = FileUpload.objects.create(
        path=source_path, name=name, user=user, import_source="S")  # ImportSource.SEQAUTO
    uploaded_vcf = UploadedVCF.objects.create(file_upload=file_upload, vcf=vcf)
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

    def test_no_sequencing_samples_falls_back_to_sheet(self):
        s = JointCalledVCFSerializer(data=self._payload("/d/trio_a.vcf.gz"))
        s.is_valid(raise_exception=True)
        joint_called_vcf = s.save()
        self.assertEqual(joint_called_vcf.sequencing_samples.count(), 0)
        self.assertEqual({ss.sample_name for ss in joint_called_vcf.get_sequencing_samples()},
                         {"S1", "S2", "S3"})
        self.assertFalse(joint_called_vcf.is_cross_run)


class CrossRunJointCalledVCFSerializerTests(TestCase):
    """A trio joint-called from samples sequenced on different runs (SACGF/variantgrid_sapath#415)."""

    def setUp(self):
        self.proband_run = _make_sequencing_run("CROSS_RUN_PROBAND")
        self.proband_sheet, _ = _make_sample_sheet(self.proband_run, ["PROBAND", "OTHER"],
                                                   sheet_hash="HSHP01")
        self.parents_run = _make_sequencing_run("CROSS_RUN_PARENTS")
        self.parents_sheet, _ = _make_sample_sheet(self.parents_run, ["MUM", "DAD"],
                                                   sheet_hash="HSHQ01")
        self.caller = VariantCaller.objects.create(name="gatk", version="4.0")

    @staticmethod
    def _member(sample_sheet, sample_name):
        return {
            "sample_sheet": {"sequencing_run": sample_sheet.sequencing_run.name,
                             "hash": sample_sheet.hash},
            "sample_name": sample_name,
        }

    def _payload(self, path, members):
        return {
            "path": path,
            "sample_sheet": {"sequencing_run": self.proband_run.name,
                             "hash": self.proband_sheet.hash},
            "variant_caller": {"name": self.caller.name, "version": self.caller.version,
                               "run_params": None},
            "sequencing_samples": members,
        }

    def _trio_members(self):
        return [self._member(self.proband_sheet, "PROBAND"),
                self._member(self.parents_sheet, "MUM"),
                self._member(self.parents_sheet, "DAD")]

    def _save(self, path, members):
        s = JointCalledVCFSerializer(data=self._payload(path, members))
        s.is_valid(raise_exception=True)
        return s.save()

    def test_members_span_two_runs(self):
        joint_called_vcf = self._save("/d/family.vcf.gz", self._trio_members())
        self.assertEqual({ss.sample_name for ss in joint_called_vcf.sequencing_samples.all()},
                         {"PROBAND", "MUM", "DAD"})
        # The owning sheet stays the one the path sits under
        self.assertEqual(joint_called_vcf.sample_sheet, self.proband_sheet)
        self.assertTrue(joint_called_vcf.is_cross_run)
        self.assertEqual(set(joint_called_vcf.sequencing_runs),
                         {self.proband_run, self.parents_run})
        # "OTHER" is on the owning sheet but is not a member
        self.assertEqual(joint_called_vcf.sequencing_sample_count, 3)

    def test_repost_replaces_members(self):
        self._save("/d/family.vcf.gz", self._trio_members())
        joint_called_vcf = self._save("/d/family.vcf.gz",
                                      [self._member(self.proband_sheet, "PROBAND"),
                                       self._member(self.parents_sheet, "MUM")])
        self.assertEqual(JointCalledVCF.objects.filter(path="/d/family.vcf.gz").count(), 1)
        self.assertEqual({ss.sample_name for ss in joint_called_vcf.sequencing_samples.all()},
                         {"PROBAND", "MUM"})

    def test_unknown_member_raises(self):
        members = self._trio_members() + [self._member(self.parents_sheet, "NOT_SEQUENCED")]
        s = JointCalledVCFSerializer(data=self._payload("/d/family.vcf.gz", members))
        s.is_valid(raise_exception=True)
        with self.assertRaises(ValidationError):
            s.save()

    def test_contributing_sheet_lists_cross_run_vcf(self):
        joint_called_vcf = self._save("/d/family.vcf.gz", self._trio_members())
        # The parents' run shows the family VCF even though it's owned by the proband's run
        self.assertIn(joint_called_vcf, self.parents_sheet.get_joint_called_vcfs())
        self.assertIn(joint_called_vcf, self.proband_sheet.get_joint_called_vcfs())


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
                                          sequencing_sample=seq_sample,
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


class CrossRunLinkSamplesTests(TestCase):
    """Linking a joint call whose members were sequenced on different runs."""

    def setUp(self):
        self.user = _make_user()
        self.caller = VariantCaller.objects.create(name="gatk", version="4.0")
        self.wgs_kit = EnrichmentKit.objects.create(name="WGS", variant_zygosity_count=False)

    def _make_run_with_member(self, run_name, sample_name, enrichment_kit=None, sheet_hash="H"):
        sequencing_run = _make_sequencing_run(run_name, enrichment_kit=enrichment_kit)
        _, sequencing_samples = _make_sample_sheet(sequencing_run, [sample_name],
                                                   sheet_hash=sheet_hash)
        return sequencing_run, sequencing_samples[0]

    def _make_trio(self, kits=(None, None, None)):
        """Three runs, one member each, joined into a single family VCF owned by the proband's run."""
        members = []
        runs = []
        for i, (sample_name, kit) in enumerate(zip(["PROBAND", "MUM", "DAD"], kits)):
            sequencing_run, member = self._make_run_with_member(
                f"XRUN_{i}", sample_name, enrichment_kit=kit, sheet_hash=f"XH{i}")
            runs.append(sequencing_run)
            members.append(member)

        path = "/d/family.vcf.gz"
        vcf, samples, uploaded_vcf = _make_vcf(self.user, "family", ["PROBAND", "MUM", "DAD"], path)
        joint_called_vcf = JointCalledVCF.objects.create(
            sequencing_run=runs[0],
            sample_sheet=members[0].sample_sheet,
            variant_caller=self.caller,
            path=path,
        )
        joint_called_vcf.sequencing_samples.set(members)
        backend = BackendVCF.objects.create(uploaded_vcf=uploaded_vcf,
                                            joint_called_vcf=joint_called_vcf,
                                            single_sample_vcf=None)
        return backend, joint_called_vcf, members, samples

    def _make_upload_step(self, backend):
        pipeline = UploadPipeline.objects.create(file_upload=backend.uploaded_vcf.file_upload)
        return UploadStep.objects.create(name="Process VCF File", upload_pipeline=pipeline,
                                         sort_order=1)

    def test_all_members_linked_across_runs(self):
        backend, _, members, samples = self._make_trio()
        link_samples_and_vcfs_to_sequencing(backend)

        linked = SampleFromSequencingSample.objects.filter(sample__in=samples)
        self.assertEqual(linked.count(), 3)
        self.assertEqual({sfss.sequencing_sample_id for sfss in linked},
                         {m.pk for m in members})

    def test_owning_run_gets_the_grid_row(self):
        backend, joint_called_vcf, _, _ = self._make_trio()
        link_samples_and_vcfs_to_sequencing(backend)

        vcf_from_run = VCFFromSequencingRun.objects.get(vcf=backend.vcf)
        self.assertEqual(vcf_from_run.sequencing_run, joint_called_vcf.sample_sheet.sequencing_run)

    def test_shared_enrichment_kit_is_applied(self):
        backend, joint_called_vcf, _, _ = self._make_trio(kits=(self.wgs_kit,) * 3)
        self.assertEqual(joint_called_vcf.enrichment_kit, self.wgs_kit)
        self.assertFalse(joint_called_vcf.has_mixed_enrichment_kits)

        link_samples_and_vcfs_to_sequencing(backend)
        backend.vcf.refresh_from_db()
        self.assertFalse(backend.vcf.variant_zygosity_count)

    def test_mixed_enrichment_kits_leave_kit_unset_and_warn(self):
        exome_kit = EnrichmentKit.objects.create(name="Exome", variant_zygosity_count=False)
        backend, joint_called_vcf, _, _ = self._make_trio(
            kits=(self.wgs_kit, exome_kit, exome_kit))
        self.assertIsNone(joint_called_vcf.enrichment_kit)
        self.assertTrue(joint_called_vcf.has_mixed_enrichment_kits)

        upload_step = self._make_upload_step(backend)
        link_samples_and_vcfs_to_sequencing(backend, upload_step=upload_step)

        backend.vcf.refresh_from_db()
        # Kit settings were skipped, so variant_zygosity_count keeps its model default
        self.assertTrue(backend.vcf.variant_zygosity_count)
        self.assertTrue(SimpleVCFImportInfo.objects.filter(
            upload_step=upload_step, message_string=MIXED_ENRICHMENT_KIT_MESSAGE).exists())

    def test_pages_show_cross_run_before_the_vcf_is_imported(self):
        """Both contributing runs surface the family VCF, and the marker shows while sample_count is None."""
        user = User.objects.create_superuser(username="crossrunadmin", password="x")
        self.client.force_login(user)

        _, joint_called_vcf, members, _ = self._make_trio()

        response = self.client.get(reverse("view_joint_called_vcf",
                                           kwargs={"joint_called_vcf_id": joint_called_vcf.pk}))
        self.assertEqual(response.status_code, 200)
        html = response.content.decode()
        self.assertIn("cross-run", html)
        for member in members:
            self.assertIn(str(member.sample_sheet.sequencing_run), html)

        for member in members:
            sequencing_run = member.sample_sheet.sequencing_run
            response = self.client.get(reverse("view_sequencing_run",
                                               kwargs={"sequencing_run_id": sequencing_run.pk}))
            self.assertEqual(response.status_code, 200, sequencing_run)
            self.assertIn("[cross-run]", response.content.decode(), sequencing_run)

    def test_new_sample_sheet_on_member_run_remaps_that_member(self):
        backend, joint_called_vcf, members, samples = self._make_trio()
        link_samples_and_vcfs_to_sequencing(backend)

        mum = members[1]
        mum_run = mum.sample_sheet.sequencing_run
        new_sheet, new_sequencing_samples = _make_sample_sheet(
            mum_run, ["MUM"], sheet_hash="XH1_v2", set_current=False)
        current = mum_run.sequencingruncurrentsamplesheet
        current_sample_sheet_changed(current, new_sheet)

        joint_called_vcf.refresh_from_db()
        member_pks = set(joint_called_vcf.sequencing_samples.values_list("pk", flat=True))
        self.assertIn(new_sequencing_samples[0].pk, member_pks)
        self.assertNotIn(mum.pk, member_pks)
        # The other two members are untouched
        self.assertIn(members[0].pk, member_pks)
        self.assertIn(members[2].pk, member_pks)
