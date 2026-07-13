from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from snpdb.models.models_enums import ImportSource
from upload.models import UploadedFile, UploadedVCF, UploadPipeline, UploadedFileTypes, ProcessingStatus
from upload.vcf.vcf_import import create_backend_vcf_links


class BackendVCFLinksTest(TestCase):
    """ 'path' is a SeqAuto-only backend-link hint - a general (non-SeqAuto) upload that supplies a
        path must not attempt (and fail) backend VCF linking. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="backend_links_user", password="x")

    def _make_uploaded_vcf(self, path):
        uploaded_file = UploadedFile.objects.create(
            user=self.user,
            name="api_upload.vcf",
            path=path,
            file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.WEB_UPLOAD,
        )
        pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, uploaded_file=uploaded_file)
        return UploadedVCF.objects.create(uploaded_file=uploaded_file, upload_pipeline=pipeline)

    @override_settings(SEQAUTO_ENABLED=False)
    def test_unmatched_path_skips_linking_when_seqauto_disabled(self):
        uploaded_vcf = self._make_uploaded_vcf(path="/client/machine/foo.vcf.gz")
        self.assertIsNone(create_backend_vcf_links(uploaded_vcf))

    @override_settings(SEQAUTO_ENABLED=True)
    def test_unmatched_path_errors_when_seqauto_enabled(self):
        uploaded_vcf = self._make_uploaded_vcf(path="/client/machine/foo.vcf.gz")
        with self.assertRaises(ValueError):
            create_backend_vcf_links(uploaded_vcf)
