import operator
from abc import abstractmethod
from functools import reduce

from library.utils import full_class_name, is_not_none
from upload.import_task_factories.import_task_factory import ImportTaskFactory
from upload.models import UploadStep, UploadStepTaskType, VCFPipelineStage
from upload.tasks.vcf.import_vcf_step_task import pipeline_start_task
from upload.tasks.vcf.import_vcf_tasks import CheckStartAnnotationTask, ScheduleMultiFileOutputTasksTask, \
    PreprocessVCFTask, UploadPipelineFinishedTask, DoNothingVCFTask


class AbstractVCFImportTaskFactory(ImportTaskFactory):

    def __init__(self):
        self.sort_order = 0

    def get_sort_order(self):
        self.sort_order += 1
        return self.sort_order

    def get_possible_extensions(self):
        return ['vcf', 'vcf.gz']

    @abstractmethod
    def get_create_data_from_vcf_header_task_class(self):
        pass

    @abstractmethod
    def get_known_variants_parallel_vcf_processing_task_class(self):
        """ Scheduled via ScheduleMultiFileOutputTasksTask after unknown variants inserted.
            One copy run for every UploadStepMultiFileOutput """
        pass

    def get_post_data_insertion_classes(self):
        return []

    def get_finish_task_classes(self):
        return [UploadPipelineFinishedTask]

    def get_pre_vcf_task(self, upload_pipeline):
        """ Run before loading initial VCF (use for eg retrieving/making it) """
        return None

    def get_create_data_from_vcf_header_task(self, upload_pipeline, input_filename):
        create_data_from_vcf_header_task_clazz = self.get_create_data_from_vcf_header_task_class()
        class_name = full_class_name(create_data_from_vcf_header_task_clazz)
        kwargs = {}
        if input_filename:
            kwargs = {"input_filename": input_filename}

        unknown_variants_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                          name="Create Data from VCF Header",
                                                          sort_order=self.get_sort_order(),
                                                          task_type=UploadStepTaskType.CELERY,
                                                          script=class_name,
                                                          **kwargs)
        return create_data_from_vcf_header_task_clazz.si(unknown_variants_step.pk, 0)

    def get_preprocess_vcf_step_and_task(self, upload_pipeline, input_filename):
        preprocess_clazz = PreprocessVCFTask
        class_name = full_class_name(preprocess_clazz)
        preprocess_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                    name=UploadStep.PREPROCESS_VCF_NAME,
                                                    sort_order=self.get_sort_order(),
                                                    task_type=UploadStepTaskType.CELERY,
                                                    pipeline_stage=VCFPipelineStage.INSERT_UNKNOWN_VARIANTS,
                                                    input_filename=input_filename,
                                                    script=class_name)
        preprocess_task = preprocess_clazz.si(preprocess_step.pk, 0)
        return preprocess_step, preprocess_task

    def create_insert_unknown_dependent_upload_steps(self, upload_pipeline, unknown_variants_step):
        """ CheckAnnotationTask decides whether to launch annotation.

            ScheduleMultiFileOutputTasksTask will launch a known_variants_parallel_vcf_processing_task_class
            for each sub file created by unknown_variants_step """
        self.create_pipeline_stage_dependent_upload_steps(upload_pipeline, VCFPipelineStage.INSERT_UNKNOWN_VARIANTS, [CheckStartAnnotationTask])

        known_variants_parallel_vcf_processing_task_class = self.get_known_variants_parallel_vcf_processing_task_class()
        schedule_class_name = full_class_name(ScheduleMultiFileOutputTasksTask)
        known_variants_parallel_vcf_processing_task_class_name = full_class_name(known_variants_parallel_vcf_processing_task_class)

        UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                  name="Schedule Parallel VCF Processing tasks",
                                  sort_order=self.get_sort_order(),
                                  task_type=UploadStepTaskType.CELERY,
                                  pipeline_stage_dependency=VCFPipelineStage.INSERT_UNKNOWN_VARIANTS,
                                  input_upload_step=unknown_variants_step,
                                  script=schedule_class_name,
                                  child_script=known_variants_parallel_vcf_processing_task_class_name)

    def create_pipeline_stage_dependent_upload_steps(self, upload_pipeline, pipeline_stage_dependency, task_classes):
        for task_class in task_classes:
            task_class_name = full_class_name(task_class)

            UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                      name=task_class.__name__,
                                      sort_order=self.get_sort_order(),
                                      task_type=UploadStepTaskType.CELERY,
                                      pipeline_stage_dependency=pipeline_stage_dependency,
                                      script=task_class_name)

    def get_fake_data_insertion_task(self, upload_pipeline):
        """ This is only needed if no other processing file tasks are run to kick off dependent steps"""
        clazz = DoNothingVCFTask
        class_name = full_class_name(clazz)
        preprocess_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                    name='Schedule tasks',
                                                    sort_order=self.get_sort_order(),
                                                    task_type=UploadStepTaskType.CELERY,
                                                    pipeline_stage=VCFPipelineStage.DATA_INSERTION,
                                                    script=class_name)
        return clazz.si(preprocess_step.pk, 0)

    def _get_vcf_filename(self, upload_pipeline) -> str:
        return upload_pipeline.uploaded_file.get_filename()

    def create_import_task(self, upload_pipeline):
        vcf_filename = self._get_vcf_filename(upload_pipeline)

        start_task = pipeline_start_task.si(upload_pipeline.pk)  # @UndefinedVariable
        pre_vcf_task = self.get_pre_vcf_task(upload_pipeline)
        create_data_from_vcf_header_task = self.get_create_data_from_vcf_header_task(upload_pipeline, vcf_filename)

        pipeline_tasks = [start_task, pre_vcf_task, create_data_from_vcf_header_task]
        if vcf_filename:
            preprocess_vcf_step, preprocess_vcf_task = self.get_preprocess_vcf_step_and_task(upload_pipeline, vcf_filename)
            pipeline_tasks.append(preprocess_vcf_task)
            # Create steps that will run when pipeline_stage is completed
            # @see upload.tasks.vcf.import_vcf_step_task.ImportVCFStepTask.check_pipeline_stage
            self.create_insert_unknown_dependent_upload_steps(upload_pipeline, preprocess_vcf_step)
            post_data_insertion_classes = self.get_post_data_insertion_classes()
        else:
            # Need one to kick off data_insertion dependencies
            pipeline_tasks.append(self.get_fake_data_insertion_task(upload_pipeline))
            post_data_insertion_classes = []

        DEPENDENT_STAGES = {
            VCFPipelineStage.DATA_INSERTION: post_data_insertion_classes,
            VCFPipelineStage.FINISH: self.get_finish_task_classes(),
        }
        for pipeline_stage_dependency, task_classes in DEPENDENT_STAGES.items():
            self.create_pipeline_stage_dependent_upload_steps(upload_pipeline, pipeline_stage_dependency, task_classes)

        pipeline_tasks = filter(is_not_none, pipeline_tasks)
        chained_pipeline_task = reduce(operator.or_, pipeline_tasks)
        return chained_pipeline_task
