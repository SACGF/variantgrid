import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.models.models import ManualVariantEntryCollection
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedFileTypes,
    UploadedManualVariantEntryCollection,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
    UploadStepTaskType,
)
from upload.tasks.vcf.import_vcf_tasks import ImportCreateUploadedVCFTask
from upload.uploaded_file_type import retry_upload_pipeline


class RetryImportTest(TestCase):
    """ Retrying an import re-runs the pipeline against records made by the original run (UploadData
        the pipeline can't re-create is kept) - so the steps have to be re-runnable """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="retry_import_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _make_pipeline(self) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.user, name="manual_variant_entry",
                                                path="/tmp/manual_variant_entry",
                                                file_type=UploadedFileTypes.MANUAL_VARIANT_ENTRY,
                                                import_source=ImportSource.WEB)
        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        UploadedManualVariantEntryCollection.objects.create(file_upload=file_upload, collection=mvec)
        return UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)

    def _make_upload_step(self, upload_pipeline: UploadPipeline) -> UploadStep:
        return UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                         name="Create Data from VCF Header",
                                         sort_order=1,
                                         task_type=UploadStepTaskType.CELERY,
                                         script="upload.tasks.vcf.import_vcf_tasks.ImportCreateUploadedVCFTask")

    def test_create_uploaded_vcf_reuses_existing(self):
        upload_step = self._make_upload_step(self._make_pipeline())
        ImportCreateUploadedVCFTask.process_items(upload_step)
        original = UploadedVCF.objects.get(upload_pipeline=upload_step.upload_pipeline)

        ImportCreateUploadedVCFTask.process_items(upload_step)  # Retry import
        self.assertEqual(UploadedVCF.objects.filter(file_upload=upload_step.file_upload).count(), 1)
        self.assertEqual(UploadedVCF.objects.get(upload_pipeline=upload_step.upload_pipeline).pk, original.pk)

    def test_retry_keeps_manual_variant_entry_collection(self):
        upload_pipeline = self._make_pipeline()
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                with patch("upload.uploaded_file_type.process_upload_pipeline") as mock_process_upload_pipeline:
                    mock_process_upload_pipeline.return_value = (upload_pipeline,)
                    retry_upload_pipeline(upload_pipeline)

        mock_process_upload_pipeline.assert_called_once()
        self.assertTrue(UploadedManualVariantEntryCollection.objects.filter(
            file_upload=upload_pipeline.file_upload).exists())
