from django.conf import settings
from django.contrib.auth.models import Group, User
from django.test import TestCase, override_settings

from seqauto.enums import DataGeneration
from seqauto.models.models_seqauto import JointCalledVCF, SampleSheet, SequencingRun
from seqauto.models.models_sequencing import Sequencer, SequencerModel
from seqauto.models.models_software import VariantCaller
from seqauto.models.models_upload import BackendVCF
from seqauto.upload_links import create_backend_vcf_links
from snpdb.enums import ImportSource
from upload.enums import UploadedFileTypes
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedVCF,
    UploadPipeline,
)
from upload.models.models import UploadUserError


class BackendVCFLinksTest(TestCase):
    """ 'path' is a SeqAuto-only backend-link hint - a general (non-SeqAuto) upload that supplies a
        path must not attempt (and fail) backend VCF linking. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="backend_links_user", password="x")

    def _make_uploaded_vcf(self, path):
        file_upload = FileUpload.objects.create(
            user=self.user,
            name="api_upload.vcf",
            path=path,
            file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.WEB_UPLOAD,
        )
        pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        return UploadedVCF.objects.create(file_upload=file_upload, upload_pipeline=pipeline)

    @override_settings(SEQAUTO_ENABLED=False)
    def test_unmatched_path_skips_linking_when_seqauto_disabled(self):
        uploaded_vcf = self._make_uploaded_vcf(path="/client/machine/foo.vcf.gz")
        self.assertIsNone(create_backend_vcf_links(uploaded_vcf))


@override_settings(SEQAUTO_ENABLED=True)
class BackendVCFLinkWriteAccessTest(TestCase):
    """ Linking an uploaded VCF to a sequencing VCF by path needs sequencing data write access, and a path
        matching no sequencing VCF fails the upload """
    PATH = "/data/RUN_2456/joint_called.vcf.gz"

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer = Sequencer.objects.create(name="MS_2456", sequencer_model=sequencer_model)
        sequencing_run = SequencingRun.objects.create(name="RUN_2456", sequencer=sequencer)
        sample_sheet = SampleSheet.objects.create(sequencing_run=sequencing_run, path="/data/RUN_2456/SampleSheet.csv")
        cls.joint_called_vcf = JointCalledVCF.objects.create(sequencing_run=sequencing_run, sample_sheet=sample_sheet,
                                                             variant_caller=VariantCaller.objects.create(name="gatk", version="4.0"),
                                                             path=cls.PATH)

    def _make_uploaded_vcf(self, user, path=PATH):
        file_upload = FileUpload.objects.create(user=user, name="api_upload.vcf", path=path,
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.WEB_UPLOAD)
        pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        return UploadedVCF.objects.create(file_upload=file_upload, upload_pipeline=pipeline)

    def test_link_refused_without_write_access(self):
        user = User.objects.create_user(username="backend_links_reader", password="x")
        uploaded_vcf = self._make_uploaded_vcf(user)
        with self.assertRaisesMessage(UploadUserError, "RUN_2456"):
            create_backend_vcf_links(uploaded_vcf)
        self.assertFalse(BackendVCF.objects.exists())

    def test_unmatched_path_fails(self):
        user = User.objects.create_user(username="backend_links_unmatched", password="x")
        user.groups.add(Group.objects.get_or_create(name=settings.SEQAUTO_API_WRITE_GROUP)[0])
        uploaded_vcf = self._make_uploaded_vcf(user, path="/client/machine/foo.vcf.gz")
        with self.assertRaisesMessage(UploadUserError, "/client/machine/foo.vcf.gz"):
            create_backend_vcf_links(uploaded_vcf)

    def test_link_made_for_sequencing_data_writer(self):
        user = User.objects.create_user(username="backend_links_writer", password="x")
        user.groups.add(Group.objects.get_or_create(name=settings.SEQAUTO_API_WRITE_GROUP)[0])
        uploaded_vcf = self._make_uploaded_vcf(user)
        backend_vcf = create_backend_vcf_links(uploaded_vcf)
        self.assertEqual(backend_vcf.joint_called_vcf, self.joint_called_vcf)
