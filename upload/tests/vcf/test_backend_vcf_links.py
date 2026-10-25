from django.conf import settings
from django.contrib.auth.models import Group, User
from django.test import TestCase, override_settings

from snpdb.enums import ImportSource
from upload.enums import UploadedFileTypes
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedVCF,
    UploadPipeline,
)
from upload.models.models import UploadUserError
from upload.vcf.vcf_import import create_backend_vcf_links


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
    def test_unmatched_path_errors_when_seqauto_enabled(self):
        self.user.groups.add(Group.objects.get_or_create(name=settings.SEQAUTO_API_WRITE_GROUP)[0])
        uploaded_vcf = self._make_uploaded_vcf(path="/client/machine/foo.vcf.gz")
        with self.assertRaises(ValueError):
            create_backend_vcf_links(uploaded_vcf)

    @override_settings(SEQAUTO_ENABLED=True)
    def test_path_needs_sequencing_data_write_access(self):
        uploaded_vcf = self._make_uploaded_vcf(path="/client/machine/foo.vcf.gz")
        with self.assertRaises(UploadUserError):
            create_backend_vcf_links(uploaded_vcf)
