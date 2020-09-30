import celery
import logging
import os

from analysis.tasks.mutational_signatures_task import calculate_mutational_signature
from annotation.annotation_versions import get_lowest_unannotated_variant_id
from annotation.models.models import AnnotationVersion
from annotation.tasks.calculate_sample_stats import calculate_vcf_stats
from library.log_utils import log_traceback
from library.utils import full_class_name
from seqauto.signals import backend_vcf_import_success_signal
from snpdb.import_status import set_vcf_and_samples_import_status
from snpdb.models import VCF
from snpdb.models.models_enums import ImportStatus, VariantsType, ProcessingStatus
from snpdb.variant_zygosity_count import update_all_variant_zygosity_counts_for_vcf, \
    create_variant_zygosity_counts
from snpdb.tasks.sample_locus_count_task import do_sample_locus_count_for_vcf_id
from upload.models import VCFPipelineStage, UploadStep, UploadStepTaskType, UploadedVCFPendingAnnotation, UploadPipeline
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.upload_processing import process_upload_pipeline
from upload.vcf.vcf_import import create_vcf_from_vcf, import_vcf_file, \
    create_import_success_message, create_cohort_genotype_collection_from_vcf
from variantgrid.celery import app


class ImportCreateVCFModelForGenotypeVCFTask(ImportVCFStepTask):
    """ Create VCF model from header """

    def process_items(self, upload_step):
        vcf_filename = upload_step.input_filename
        upload_pipeline = upload_step.upload_pipeline

        try:
            vcf = upload_pipeline.uploadedvcf.vcf
            _ = vcf.pk  # Throws exception if VCF is None
            logging.info("VCF already existed - reusing!")
        except AttributeError:
            vcf = create_vcf_from_vcf(upload_step, vcf_filename)

        # If build not set, end
        if vcf.genome_build is None:
            vcf.import_status = ImportStatus.REQUIRES_USER_INPUT
            vcf.save()

            # Make other tasks get skipped...
            upload_pipeline.status = ProcessingStatus.TERMINATED_EARLY
            upload_pipeline.save()
            return

        # If past this point - we have VCF w/genome build
        create_cohort_genotype_collection_from_vcf(vcf)

        # CalculateVCFStatsTask - called after annotation finished
        vcf_stats_task_class_name = full_class_name(CalculateVCFStatsTask)
        calculate_mutational_signature_task_class_name = full_class_name(CalculateMutationalSignaturesForSampleTask)

        sort_order = upload_step.upload_pipeline.get_max_step_sort_order()
        sort_order += 1
        UploadStep.objects.create(upload_pipeline=upload_step.upload_pipeline,
                                  name="Calculate VCF Stats",
                                  sort_order=sort_order,
                                  task_type=UploadStepTaskType.CELERY,
                                  pipeline_stage_dependency=VCFPipelineStage.ANNOTATION_COMPLETE,
                                  script=vcf_stats_task_class_name)

        for sample in vcf.sample_set.all():
            if sample.variants_type == VariantsType.SOMATIC_ONLY:
                sort_order += 1
                UploadStep.objects.create(upload_pipeline=upload_step.upload_pipeline,
                                          name="Calculate Mutational Signature",
                                          sort_order=sort_order,
                                          task_type=UploadStepTaskType.CELERY,
                                          pipeline_stage_dependency=VCFPipelineStage.ANNOTATION_COMPLETE,
                                          input_filename=str(sample.pk),
                                          script=calculate_mutational_signature_task_class_name)


class VCFCheckAnnotationTask(ImportVCFStepTask):
    """ Check if all of our variants are annotated, and trigger annotation stage steps if so
        Otherwise, will leave a UploadedVCFPendingAnnotation object with finished=NULL
        which will be tested every time 'annotation_run_complete_signal' fires """

    def process_items(self, upload_step):
        uploaded_vcf = upload_step.get_uploaded_vcf()
        if uploaded_vcf.max_variant is None:
            msg = "UploadedVCF.max_variant is None! This must be set so we know if we need to perform annotation"
            raise ValueError(msg)

        pending_annotation = UploadedVCFPendingAnnotation.objects.create(uploaded_vcf=uploaded_vcf)
        annotation_version = AnnotationVersion.latest(upload_step.genome_build)
        variant_annotation_version = annotation_version.variant_annotation_version
        lowest_unannotated_variant_id = get_lowest_unannotated_variant_id(variant_annotation_version)
        pending_annotation.attempt_schedule_annotation_stage_steps(lowest_unannotated_variant_id)


class ProcessGenotypeVCFDataTask(ImportVCFStepTask):
    """ Processes a VCF for which vcf data from header has already been made
        (ie via ImportGenotypeVCFTask) - this can run in parallel """

    def process_items(self, upload_step):
        items_processed = import_vcf_file(upload_step)
        return items_processed


class CalculateVCFStatsTask(ImportVCFStepTask):
    """ Called when pipeline_stage ANNOTATION_COMPLETE done """

    def process_items(self, upload_step):
        annotation_version = AnnotationVersion.latest(upload_step.genome_build)
        vcf = upload_step.upload_pipeline.uploadedvcf.vcf
        calculate_vcf_stats(vcf.pk, annotation_version.pk)


class CalculateMutationalSignaturesForSampleTask(ImportVCFStepTask):
    """ Called when pipeline_stage ANNOTATION_COMPLETE done """

    def process_items(self, upload_step):
        sample_id = upload_step.input_filename
        calculate_mutational_signature(sample_id)


class UpdateVariantZygosityCountsTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        uploaded_vcf = upload_step.get_uploaded_vcf()
        create_variant_zygosity_counts()
        update_all_variant_zygosity_counts_for_vcf(uploaded_vcf.vcf, '+')


class SampleLocusCountsTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        uploaded_vcf = upload_step.get_uploaded_vcf()

        do_sample_locus_count_for_vcf_id(uploaded_vcf.vcf.pk)


class ImportGenotypeVCFSuccessTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        uploaded_vcf = upload_step.get_uploaded_vcf()
        vcf = uploaded_vcf.vcf
        print(f"ImportVCFSuccessTask for VCF = {vcf}")

        set_vcf_and_samples_import_status(vcf, ImportStatus.SUCCESS)

        try:
            backend_vcf = vcf.uploadedvcf.backendvcf
            backend_vcf_import_success_signal.send(sender=os.path.basename(__file__), backend_vcf=backend_vcf)
        except:
            pass

        create_import_success_message(vcf)


@celery.task
def reload_vcf_task(upload_pipeline_id, vcf_id):
    """ Needs to be done thread safely for zygosity count """

    try:
        upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)
        if vcf_id:
            vcf = VCF.objects.get(pk=vcf_id)
            vcf.import_status = ImportStatus.IMPORTING
            vcf.save()

            logging.info("VariantCount")
            try:
                update_all_variant_zygosity_counts_for_vcf(vcf, '-')
            except:
                pass

            logging.info("Delete internal data")
            vcf.delete_internal_data()

        logging.info("Reloading pipeline")
        process_upload_pipeline(upload_pipeline)
    except:
        log_traceback()


ImportCreateVCFModelForGenotypeVCFTask = app.register_task(ImportCreateVCFModelForGenotypeVCFTask())
VCFCheckAnnotationTask = app.register_task(VCFCheckAnnotationTask())
ProcessGenotypeVCFDataTask = app.register_task(ProcessGenotypeVCFDataTask())
CalculateVCFStatsTask = app.register_task(CalculateVCFStatsTask())
CalculateMutationalSignaturesForSampleTask = app.register_task(CalculateMutationalSignaturesForSampleTask())
UpdateVariantZygosityCountsTask = app.register_task(UpdateVariantZygosityCountsTask())
SampleLocusCountsTask = app.register_task(SampleLocusCountsTask())
ImportGenotypeVCFSuccessTask = app.register_task(ImportGenotypeVCFSuccessTask())
