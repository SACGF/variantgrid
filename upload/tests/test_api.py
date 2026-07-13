import gzip
import os
import tempfile
from hashlib import sha256
from unittest.mock import patch

from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from analysis.models import AnalysisTemplate
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import CachedGeneratedFile, GenomeBuild
from snpdb.models.models_enums import ImportSource, ProcessingStatus
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from upload.models import UploadedFile, UploadedVCF, UploadPipeline, UploadedFileTypes

COHORT_EXPORT_TEMPLATE_NAME = "Cohort VCF Export auto analysis"  # settings.ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT


def _vcf_bytes(marker: str) -> bytes:
    return (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t123\t{marker}\tA\tG\t.\t.\t.\n"
    ).encode()


class UploadAPITestBase(APITestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="api_owner", password="x")
        cls.other_user = User.objects.create_user(username="api_other", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls._temp_upload_paths = []

    @classmethod
    def tearDownClass(cls):
        for path in getattr(cls, "_temp_upload_paths", []):
            try:
                if os.path.exists(path):
                    os.unlink(path)
            except OSError:
                pass
        super().tearDownClass()

    @classmethod
    def _create_vcf_upload(cls, marker: str, *, pipeline_status=ProcessingStatus.SUCCESS, path=None):
        """ Build a fully-processed VCF upload (UploadedFile + UploadPipeline + UploadedVCF)
            attached to a permissioned cohort's VCF. """
        cohort = create_fake_cohort(cls.owner, cls.grch37)
        vcf = cohort.vcf
        django_file = SimpleUploadedFile(f"{marker}.vcf", _vcf_bytes(marker))
        uploaded_file = UploadedFile.objects.create(
            user=cls.owner,
            name=f"{marker}.vcf",
            path=path,
            file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.WEB_UPLOAD,
            uploaded_file=django_file,
        )
        uploaded_file.store_sha256_hash()
        cls._temp_upload_paths.append(uploaded_file.get_filename())
        pipeline = UploadPipeline.objects.create(status=pipeline_status, uploaded_file=uploaded_file)
        UploadedVCF.objects.create(uploaded_file=uploaded_file, vcf=vcf, upload_pipeline=pipeline)
        return uploaded_file


class UploadFileAPITest(UploadAPITestBase):
    def test_upload_populates_sha256(self):
        """ handle_file_upload stores sha256_hash for new uploads (dedup/poll by content) """
        self.client.force_authenticate(user=self.owner)
        content = _vcf_bytes("newupload")
        expected_hash = sha256(content).hexdigest()

        with patch("upload.views.views.upload_processing.process_uploaded_file"):
            response = self.client.post(
                reverse("api_file_upload"),
                {"file": SimpleUploadedFile("newupload.vcf", content)},
                format="multipart",
            )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        uploaded_file = UploadedFile.objects.get(pk=data["uploaded_file_id"])
        self._temp_upload_paths.append(uploaded_file.get_filename())
        self.assertEqual(uploaded_file.sha256_hash, expected_hash)
        self.assertEqual(data["sha256_hash"], expected_hash)

    def test_upload_dedups_by_hash_without_path(self):
        existing = self._create_vcf_upload("dedup_nopath")
        self.client.force_authenticate(user=self.owner)

        response = self.client.post(
            reverse("api_file_upload"),
            {"file": SimpleUploadedFile("dedup_nopath.vcf", _vcf_bytes("dedup_nopath"))},
            format="multipart",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertEqual(data["uploaded_file_id"], existing.pk)
        self.assertIn("message", data)

    def test_upload_dedups_by_hash_with_path(self):
        existing = self._create_vcf_upload("dedup_path", path="/client/sample.vcf")
        self.client.force_authenticate(user=self.owner)

        response = self.client.post(
            reverse("api_file_upload") + "?path=/client/sample.vcf",
            {"file": SimpleUploadedFile("dedup_path.vcf", _vcf_bytes("dedup_path"))},
            format="multipart",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()["uploaded_file_id"], existing.pk)


class UploadStatusAPITest(UploadAPITestBase):
    def test_status_by_id_processing(self):
        uploaded_file = self._create_vcf_upload("status_proc", pipeline_status=ProcessingStatus.PROCESSING)
        self.client.force_authenticate(user=self.owner)

        response = self.client.get(reverse("api_upload_status",
                                           kwargs={"uploaded_file_id": uploaded_file.pk}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertEqual(data["uploaded_file_id"], uploaded_file.pk)
        self.assertEqual(data["pipeline_status"], ProcessingStatus.PROCESSING.label)
        self.assertFalse(data["annotation_complete"])
        self.assertIsNotNone(data["vcf_id"])
        self.assertTrue(data["samples"])

    def test_status_before_vcf_created(self):
        """ Race: client polls immediately after upload, UploadedVCF exists but its vcf is still None """
        django_file = SimpleUploadedFile("racey.vcf", _vcf_bytes("racey"))
        uploaded_file = UploadedFile.objects.create(
            user=self.owner, name="racey.vcf", file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.WEB_UPLOAD, uploaded_file=django_file)
        uploaded_file.store_sha256_hash()
        self._temp_upload_paths.append(uploaded_file.get_filename())
        pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, uploaded_file=uploaded_file)
        UploadedVCF.objects.create(uploaded_file=uploaded_file, vcf=None, upload_pipeline=pipeline)

        self.client.force_authenticate(user=self.owner)
        response = self.client.get(reverse("api_upload_status",
                                           kwargs={"uploaded_file_id": uploaded_file.pk}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertIsNone(data["vcf_id"])
        self.assertEqual(data["samples"], [])
        self.assertFalse(data["annotation_complete"])

    def test_status_by_sha256_success(self):
        uploaded_file = self._create_vcf_upload("status_ok", pipeline_status=ProcessingStatus.SUCCESS)
        self.client.force_authenticate(user=self.owner)

        response = self.client.get(reverse("api_upload_status_sha256",
                                           kwargs={"sha256_hash": uploaded_file.sha256_hash}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertEqual(data["uploaded_file_id"], uploaded_file.pk)
        self.assertEqual(data["remaining_annotation_runs"], 0)
        self.assertTrue(data["annotation_complete"])
        # Default settings name a template that doesn't exist in tests -> downloads not available
        self.assertFalse(data["downloads_available"])

    def test_status_permission_denied_for_other_user(self):
        uploaded_file = self._create_vcf_upload("status_perm")
        self.client.force_authenticate(user=self.other_user)

        response = self.client.get(reverse("api_upload_status",
                                           kwargs={"uploaded_file_id": uploaded_file.pk}))
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)


class AnnotatedDownloadAPITest(UploadAPITestBase):
    def _make_export_template(self):
        return AnalysisTemplate.objects.create(name=COHORT_EXPORT_TEMPLATE_NAME, user=self.owner)

    def test_invalid_export_type(self):
        uploaded_file = self._create_vcf_upload("dl_badtype")
        self.client.force_authenticate(user=self.owner)
        response = self.client.get(reverse("api_annotated_download",
                                           kwargs={"uploaded_file_id": uploaded_file.pk,
                                                   "export_type": "bam"}))
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    def test_download_template_not_configured(self):
        uploaded_file = self._create_vcf_upload("dl_notemplate")
        self.client.force_authenticate(user=self.owner)
        response = self.client.get(reverse("api_annotated_download",
                                           kwargs={"uploaded_file_id": uploaded_file.pk,
                                                   "export_type": "vcf"}))
        self.assertEqual(response.status_code, status.HTTP_501_NOT_IMPLEMENTED)
        self.assertIn("error", response.json())

    def test_download_generating_returns_202(self):
        uploaded_file = self._create_vcf_upload("dl_generating")
        self._make_export_template()
        self.client.force_authenticate(user=self.owner)

        cgf = CachedGeneratedFile(generator="export_cohort_to_downloadable_file",
                                  params_hash="fake", progress=0.3, filename=None)
        with patch.object(CachedGeneratedFile, "get_or_create_and_launch", return_value=cgf):
            response = self.client.get(reverse("api_annotated_download",
                                               kwargs={"uploaded_file_id": uploaded_file.pk,
                                                       "export_type": "vcf"}))
        self.assertEqual(response.status_code, status.HTTP_202_ACCEPTED)
        data = response.json()
        self.assertEqual(data["status"], "generating")
        self.assertEqual(data["progress"], 0.3)

    def test_download_ready_streams_file(self):
        uploaded_file = self._create_vcf_upload("dl_ready")
        self._make_export_template()
        self.client.force_authenticate(user=self.owner)

        payload = _vcf_bytes("annotated")
        fd, gz_path = tempfile.mkstemp(suffix=".vcf.gz")
        os.close(fd)
        self._temp_upload_paths.append(gz_path)
        with gzip.open(gz_path, "wb") as f:
            f.write(payload)

        cgf = CachedGeneratedFile(generator="export_cohort_to_downloadable_file",
                                  params_hash="fake", progress=1.0, filename=gz_path)
        with patch.object(CachedGeneratedFile, "get_or_create_and_launch", return_value=cgf):
            response = self.client.get(reverse("api_annotated_download",
                                               kwargs={"uploaded_file_id": uploaded_file.pk,
                                                       "export_type": "vcf"}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response["Content-Type"], "application/gzip")
        self.assertIn("attachment", response["Content-Disposition"])
        streamed = b"".join(response.streaming_content)
        self.assertEqual(gzip.decompress(streamed), payload)
