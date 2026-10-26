"""
FINISH is scheduled once per DATA_INSERTION step that sees it is last, so it can run more than once - only
the call that claims the FINISH steps may queue the chain that ends in pipeline_success_task.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from library.utils import full_class_name
from snpdb.enums import ImportSource
from upload.enums import UploadedFileTypes, UploadStepTaskType, VCFPipelineStage
from upload.models import FileUpload, ProcessingStatus, UploadPipeline, UploadStep
from upload.tasks.vcf.genotype_vcf_tasks import ImportGenotypeVCFSuccessTask, reload_vcf_task
from upload.tasks.vcf.import_vcf_step_task import schedule_pipeline_stage_steps


class SchedulePipelineStageStepsTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="schedule_stage_user", password="x")

    def _make_pipeline(self, status=ProcessingStatus.PROCESSING) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.user, name="x.vcf", path="/tmp/x.vcf",
                                                file_type=UploadedFileTypes.VCF, import_source=ImportSource.WEB)
        return UploadPipeline.objects.create(status=status, file_upload=file_upload)

    @patch("upload.tasks.vcf.import_vcf_step_task.chain")
    def test_repeat_finish_queues_nothing(self, mock_chain):
        upload_pipeline = self._make_pipeline()
        UploadStep.objects.create(upload_pipeline=upload_pipeline, name="Success", sort_order=1,
                                  task_type=UploadStepTaskType.CELERY,
                                  pipeline_stage_dependency=VCFPipelineStage.FINISH,
                                  script=full_class_name(ImportGenotypeVCFSuccessTask))

        schedule_pipeline_stage_steps(upload_pipeline.pk, VCFPipelineStage.FINISH)
        schedule_pipeline_stage_steps(upload_pipeline.pk, VCFPipelineStage.FINISH)
        self.assertEqual(mock_chain.call_count, 1)
        self.assertEqual(len(mock_chain.call_args.args[0]), 2)  # The step, then pipeline_success_task

    @patch("upload.tasks.vcf.import_vcf_step_task.chain")
    def test_finish_without_steps_still_closes(self, mock_chain):
        upload_pipeline = self._make_pipeline()
        schedule_pipeline_stage_steps(upload_pipeline.pk, VCFPipelineStage.FINISH)
        self.assertEqual(mock_chain.call_count, 1)

    @patch("upload.tasks.vcf.genotype_vcf_tasks.process_upload_pipeline", side_effect=ValueError("boom"))
    def test_reload_failure_leaves_pipeline_retryable(self, _process_upload_pipeline):
        upload_pipeline = self._make_pipeline(status=ProcessingStatus.CREATED)
        with patch("upload.models.models.report_exc_info"), patch("upload.models.models.report_message"):
            reload_vcf_task(upload_pipeline.pk, None)
        upload_pipeline.refresh_from_db()
        self.assertEqual(upload_pipeline.status, ProcessingStatus.ERROR)
        self.assertIn("boom", upload_pipeline.progress_status)
