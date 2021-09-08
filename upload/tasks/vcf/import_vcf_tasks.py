import logging

import celery
from django.contrib.auth.models import User

from annotation.tasks.annotation_scheduler_task import annotation_scheduler
from library.log_utils import log_traceback
from library.utils import full_class_name, import_class
from snpdb.models.models_enums import ProcessingStatus, AlleleConversionTool
from snpdb.ncbi_remap import ncbi_remap
from upload.models import UploadStep, UploadStepTaskType, VCFPipelineStage, UploadedVCF
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.upload_processing import process_vcf_file
from upload.vcf.bulk_allele_linking_vcf_processor import BulkAlleleLinkingVCFProcessor
from upload.vcf.bulk_minimal_vcf_processor import BulkMinimalVCFProcessor
from upload.vcf.vcf_import import import_vcf_file, get_preprocess_vcf_import_info
from upload.vcf.vcf_preprocess import preprocess_vcf
from variantgrid.celery import app


class PreprocessVCFTask(ImportVCFStepTask):
    """ * Decompose (multi-alts to different lines) using vt
        * Normalise (@see https://genome.sph.umich.edu/wiki/Variant_Normalization) using vt
        * Split VCF file into sub files (so it can be processed in parallel)

        This VCFStepTask has multi-file output, get via upload_step.get_input_filenames() """

    def process_items(self, upload_step):
        preprocess_vcf(upload_step)

        # Reload from DB - vcf_extract_unknown_and_split_file set items_processed in a different process
        upload_step = UploadStep.objects.get(pk=upload_step.pk)
        return upload_step.items_processed


class CheckStartAnnotationTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        if upload_step.pipeline_inserted_unknown_variants():
            annotation_scheduler.apply_async()  # @UndefinedVariable
        return 0


class ScheduleMultiFileOutputTasksTask(ImportVCFStepTask):
    """ Creates and launches tasks of type upload_step.child_script for each
        multi-file output from upload_step.input_upload_step

        Example is to split a VCF up into multiple sub files, then launch parallel processors """

    def process_items(self, upload_step):
        child_task_class = import_class(upload_step.child_script)
        child_class_name = full_class_name(child_task_class)
        sort_order = upload_step.upload_pipeline.get_max_step_sort_order()

        for input_filename, items_to_process in upload_step.get_multi_input_files_and_records():
            sort_order += 1
            child_step = UploadStep.objects.create(upload_pipeline=upload_step.upload_pipeline,
                                                   name=UploadStep.PROCESS_VCF_TASK_NAME,
                                                   sort_order=sort_order,
                                                   task_type=UploadStepTaskType.CELERY,
                                                   pipeline_stage=VCFPipelineStage.DATA_INSERTION,
                                                   script=child_class_name,
                                                   input_filename=input_filename,
                                                   items_to_process=items_to_process)
            child_step.launch_task(child_task_class)

        return 0


class UploadPipelineFinishedTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        upload_pipeline = upload_step.upload_pipeline
        if upload_pipeline.status == ProcessingStatus.PROCESSING:
            upload_pipeline.status = ProcessingStatus.SUCCESS
            upload_pipeline.save()
        return 0


class ImportCreateUploadedVCFTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        upload_pipeline = upload_step.upload_pipeline
        uploaded_file = upload_pipeline.uploaded_file
        UploadedVCF.objects.create(uploaded_file=uploaded_file,
                                   upload_pipeline=upload_pipeline)
        return 0


class ProcessVCFSetMaxVariantTask(ImportVCFStepTask):
    """ Finds highest variant_id in VCF so we can tell whether we're done annotating or not
        Can run in parallel on split VCFs """

    def process_items(self, upload_step):
        preprocess_vcf_import_info = get_preprocess_vcf_import_info(upload_step.upload_pipeline)
        bulk_inserter = BulkMinimalVCFProcessor(upload_step, preprocess_vcf_import_info)
        items_processed = import_vcf_file(upload_step, bulk_inserter)
        return items_processed


class ProcessVCFLinkAllelesSetMaxVariantTask(ImportVCFStepTask):
    """ Link Alleles provided as the ID column in VCF
        Finds highest variant_id in VCF so we can tell whether we're done annotating or not
        Can run in parallel on split VCFs """

    def process_items(self, upload_step):
        preprocess_vcf_import_info = get_preprocess_vcf_import_info(upload_step.upload_pipeline)
        bulk_inserter = BulkAlleleLinkingVCFProcessor(upload_step, preprocess_vcf_import_info)
        items_processed = import_vcf_file(upload_step, bulk_inserter)
        return items_processed


class DoNothingVCFTask(ImportVCFStepTask):
    """ Used to be able to check dependencies """

    def process_items(self, upload_step):
        pass


class LiftoverCreateVCFTask(ImportVCFStepTask):
    def process_items(self, upload_step: UploadStep):
        upload_pipeline = upload_step.upload_pipeline
        liftover = upload_pipeline.uploaded_file.uploadedliftover.liftover
        if liftover.conversion_tool == AlleleConversionTool.NCBI_REMAP:
            ncbi_remap(upload_step.input_filename, liftover.source_genome_build,
                       upload_step.output_filename, liftover.genome_build)
        else:
            raise ValueError(f"Don't know how to produce liftover VCF for {liftover.get_conversion_tool_display()}")


class LiftoverCompleteTask(ImportVCFStepTask):

    def process_items(self, upload_step: UploadStep):
        upload_pipeline = upload_step.upload_pipeline
        liftover = upload_pipeline.uploaded_file.uploadedliftover.liftover
        liftover.complete()


@celery.task
def process_vcf_file_task(vcf_filename, name, user_id, import_source):
    """ This has been made a task so we can call it by name (not imports) via seqauto """

    try:
        user = User.objects.get(pk=user_id)
        process_vcf_file(vcf_filename, name, user, import_source=import_source)
        logging.info("process_vcf_file_task: Finishing normally")
    except:
        logging.warning("process_vcf_file_task: Got exception")
        log_traceback()


PreprocessVCFTask = app.register_task(PreprocessVCFTask())
CheckStartAnnotationTask = app.register_task(CheckStartAnnotationTask())
ScheduleMultiFileOutputTasksTask = app.register_task(ScheduleMultiFileOutputTasksTask())
UploadPipelineFinishedTask = app.register_task(UploadPipelineFinishedTask())
ImportCreateUploadedVCFTask = app.register_task(ImportCreateUploadedVCFTask())
ProcessVCFSetMaxVariantTask = app.register_task(ProcessVCFSetMaxVariantTask())
ProcessVCFLinkAllelesSetMaxVariantTask = app.register_task(ProcessVCFLinkAllelesSetMaxVariantTask())
DoNothingVCFTask = app.register_task(DoNothingVCFTask())
LiftoverCreateVCFTask = app.register_task(LiftoverCreateVCFTask())
LiftoverCompleteTask = app.register_task(LiftoverCompleteTask())
