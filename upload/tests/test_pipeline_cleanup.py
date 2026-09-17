"""
#928: a VCF pipeline is closed in exactly one place, and its scratch files go with it.

Covers the two things that leaked - pipeline_success_task being defeated by a FINISH step that set
SUCCESS itself, and nothing removing a pipeline's directories when its row was deleted.
"""
import os
import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.models.models import ManualVariantEntryCollection
from library.django_utils.django_file_utils import (
    get_import_processing_dir,
    import_processing_dir_path,
)
from snpdb.clingen_allele_api import ClinGenAlleleRegistryAPI
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedFileTypes,
    UploadedManualVariantEntryCollection,
    UploadPipeline,
    UploadStep,
    UploadStepTaskType,
    VCFPipelineStage,
)
from upload.tasks.vcf.import_vcf_step_task import pipeline_success_task
from upload.uploaded_file_type import retry_upload_pipeline


class ImportProcessingDirTest(TestCase):

    def test_path_helper_creates_nothing(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                path = import_processing_dir_path(42)
                self.assertFalse(os.path.exists(path))
                self.assertTrue(os.path.exists(get_import_processing_dir(42)))

    def test_clingen_api_creates_nothing(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                ClinGenAlleleRegistryAPI()
                self.assertEqual(os.listdir(import_processing_dir), [])


class PipelineCleanupTest(TestCase):
    """ Manual variant entry is the smallest pipeline whose input VCF we generate ourselves """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="pipeline_cleanup_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _make_pipeline(self, generated_input=True) -> UploadPipeline:
        """ Writes the generated input VCF into its own manual_variants_<pk> dir, as
            annotation.manual_variant_entry.create_manual_variants does """
        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        if generated_input:
            working_dir = get_import_processing_dir(mvec.pk, "manual_variants")
            path = os.path.join(working_dir, "manual_variant_entry.vcf")
        else:
            path = os.path.join(tempfile.gettempdir(), "somebody_elses.vcf")
        with open(path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")

        file_upload = FileUpload.objects.create(user=self.user, name="manual_variant_entry", path=path,
                                                file_type=UploadedFileTypes.MANUAL_VARIANT_ENTRY,
                                                import_source=ImportSource.WEB)
        UploadedManualVariantEntryCollection.objects.create(file_upload=file_upload, collection=mvec)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING,
                                                        file_upload=file_upload)
        # Something in the pipeline's own dir, so we can see it go
        with open(os.path.join(upload_pipeline.get_pipeline_processing_dir(), "split_1.vcf"), "w") as f:
            f.write("")
        return upload_pipeline

    def _make_finish_step(self, upload_pipeline: UploadPipeline) -> UploadStep:
        return UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                         name="LiftoverCompleteTask",
                                         sort_order=1,
                                         task_type=UploadStepTaskType.CELERY,
                                         pipeline_stage_dependency=VCFPipelineStage.FINISH,
                                         script="upload.tasks.vcf.import_vcf_tasks.LiftoverCompleteTask",
                                         start_date="2026-09-08T00:00:00Z",
                                         end_date="2026-09-08T00:01:00Z")

    def test_success_closes_pipeline_and_removes_files(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                upload_pipeline = self._make_pipeline()
                self._make_finish_step(upload_pipeline)
                pipeline_dir = import_processing_dir_path(upload_pipeline.pk)

                pipeline_success_task(upload_pipeline.pk)

                upload_pipeline.refresh_from_db()
                self.assertEqual(upload_pipeline.status, ProcessingStatus.SUCCESS)
                self.assertEqual(upload_pipeline.processing_seconds_wall_time, 60)
                self.assertFalse(os.path.exists(pipeline_dir))
                self.assertFalse(os.path.exists(upload_pipeline.file_upload.path))

    def test_success_keeps_files_when_setting_off(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir,
                                   IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS=False):
                upload_pipeline = self._make_pipeline()
                self._make_finish_step(upload_pipeline)

                pipeline_success_task(upload_pipeline.pk)

                upload_pipeline.refresh_from_db()
                self.assertEqual(upload_pipeline.status, ProcessingStatus.SUCCESS)
                self.assertTrue(os.path.exists(import_processing_dir_path(upload_pipeline.pk)))
                self.assertTrue(os.path.exists(upload_pipeline.file_upload.path))

    def test_success_task_is_a_no_op_once_closed(self):
        """ FINISH can be scheduled more than once, so the second run has to leave a closed pipeline
            alone rather than re-fire its success event """
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                upload_pipeline = self._make_pipeline()
                self._make_finish_step(upload_pipeline)
                pipeline_success_task(upload_pipeline.pk)
                UploadPipeline.objects.filter(pk=upload_pipeline.pk).update(items_processed=123)

                pipeline_success_task(upload_pipeline.pk)

                upload_pipeline.refresh_from_db()
                self.assertEqual(upload_pipeline.status, ProcessingStatus.SUCCESS)
                self.assertEqual(upload_pipeline.items_processed, 123)

    def test_delete_removes_files(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir,
                                   IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS=False):
                upload_pipeline = self._make_pipeline()
                pipeline_dir = import_processing_dir_path(upload_pipeline.pk)
                input_path = upload_pipeline.file_upload.path

                upload_pipeline.delete()

                self.assertFalse(os.path.exists(pipeline_dir))
                self.assertFalse(os.path.exists(input_path))

    def test_file_upload_cascade_removes_files(self):
        """ Deleting the FileUpload (eg a user deleting their upload) cascades to the pipeline """
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                upload_pipeline = self._make_pipeline()
                pipeline_dir = import_processing_dir_path(upload_pipeline.pk)
                input_path = upload_pipeline.file_upload.path

                upload_pipeline.file_upload.delete()

                self.assertFalse(os.path.exists(pipeline_dir))
                self.assertFalse(os.path.exists(input_path))

    def test_delete_keeps_a_file_we_did_not_generate(self):
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                upload_pipeline = self._make_pipeline(generated_input=False)
                input_path = upload_pipeline.file_upload.path
                try:
                    upload_pipeline.delete()
                    self.assertTrue(os.path.exists(input_path))
                finally:
                    os.remove(input_path)

    def test_retry_keeps_the_generated_input(self):
        """ Retry re-runs the pipeline against the VCF we generated for it, so only the pipeline's own
            working files go """
        with tempfile.TemporaryDirectory() as import_processing_dir:
            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                upload_pipeline = self._make_pipeline()
                pipeline_dir = import_processing_dir_path(upload_pipeline.pk)
                input_path = upload_pipeline.file_upload.path

                with patch("upload.uploaded_file_type.process_upload_pipeline") as mock_process:
                    mock_process.return_value = (upload_pipeline,)
                    retry_upload_pipeline(upload_pipeline)

                self.assertFalse(os.path.exists(pipeline_dir))
                self.assertTrue(os.path.exists(input_path))
