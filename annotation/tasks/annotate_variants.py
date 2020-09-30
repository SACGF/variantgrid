from celery import chain
from django.conf import settings
from django.utils import timezone
import celery
import logging
import os

from annotation.models import AnnotationStatus
from annotation.models.models import AnnotationRun
from annotation.signals import annotation_run_complete_signal
from annotation.vcf_files.variants_to_vcf import unannotated_variants_to_vcf
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import get_vep_command
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.file_utils import name_from_filename
from library.log_utils import get_traceback, report_message
from library.utils import execute_cmd


@celery.task
def delete_annotation_run(annotation_run_id):
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    annotation_run.status = AnnotationStatus.DELETING
    annotation_run.save()
    annotation_run.delete()


@celery.task
def delete_annotation_run_uploaded_data(annotation_run_id):
    """ Deletes related objects but not actual run (used for retry annotation upload) """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    annotation_run.delete_related_objects()


@celery.task
def assign_range_lock_to_annotation_run(annotation_run_id, annotation_range_lock_id):
    # has 1-to-1 so any previous AnnotationRuns linking to AnnotationRangeLock must have been deleted
    AnnotationRun.objects.filter(pk=annotation_run_id).update(annotation_range_lock_id=annotation_range_lock_id)


@celery.task
def annotate_variants(annotation_run_id):
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    logging.info("annotate_variants: %s", annotation_run)

    # task_id used as Celery lock
    num_modified = AnnotationRun.objects.filter(pk=annotation_run.pk,
                                                task_id__isnull=True).update(task_id=annotate_variants.request.id)
    if num_modified != 1:
        msg = f"Celery could't get task_id lock on AnnotationRun: {annotation_run.pk}"
        annotation_run.celery_task_logs[annotate_variants.request.id] = msg
        annotation_run.save()
        raise ValueError(msg)

    try:
        # Reload to get updated task_id
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        annotation_run.task_id = annotate_variants.request.id
        annotation_run.set_task_log("start", timezone.now())
        annotation_run.save()

        if annotation_run.vcf_annotated_filename is None:
            dump_and_annotate_variants(annotation_run)

        if annotation_run.vcf_annotated_filename:
            import_vcf_annotations(annotation_run)

        annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                            variant_annotation_version=annotation_run.annotation_range_lock.version)
    except Exception as e:
        tb = get_traceback()
        error_message = f"{e}: {annotation_run.pipeline_stderr}"

        name = 'Annotation pipeline run ' + str(annotation_run.id)
        report_message(message=name,
                       level='error',
                       extra_data={'output': error_message,
                                   'error': tb})

        annotation_run.error_exception = tb
        annotation_run.set_task_log("error_exception", tb)
        annotation_run.save()
        create_event(None, "AnnotationRun failed", tb, severity=LogLevel.ERROR)
        raise
    finally:
        annotation_run.set_task_log("end", timezone.now())
        annotation_run.task_id = None
        annotation_run.save()


def dump_and_annotate_variants(annotation_run):
    vcf_base_name = f"dump_{annotation_run.pk}.vcf"
    vcf_dump_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, vcf_base_name)

    annotation_run.dump_start = timezone.now()
    annotation_run.vcf_dump_filename = vcf_dump_filename
    annotation_run.save()

    genome_build = annotation_run.genome_build
    annotation_consortium = annotation_run.annotation_consortium
    variants_to_annotate = unannotated_variants_to_vcf(genome_build, vcf_dump_filename,
                                                       annotation_run.annotation_range_lock)

    logging.info("Annotating %d variants", variants_to_annotate)
    if variants_to_annotate:
        annotation_run.dump_end = timezone.now()
        annotation_run.save()

        name = name_from_filename(vcf_dump_filename)
        vcf_annotated_basename = f"{name}.vep_annotated_{genome_build.name}.vcf.gz"
        vcf_annotated_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, vcf_annotated_basename)

        cmd = get_vep_command(vcf_dump_filename, vcf_annotated_filename, genome_build, annotation_consortium)
        annotation_run.annotation_start = timezone.now()
        annotation_run.pipeline_command = " ".join(cmd)
        annotation_run.save()
        return_code, std_out, std_err = execute_cmd(cmd)
        annotation_run.pipeline_stdout = std_out
        annotation_run.pipeline_stderr = std_err
        logging.info(f"VEP returned code: {return_code}")

        if return_code != 0:
            annotation_run.save()  # save stdout/stderr
            raise RuntimeError(f"VEP returned {return_code}")

        now = timezone.now()
        annotation_run.vcf_annotated_filename = vcf_annotated_filename
        annotation_run.annotation_end = now
        annotation_run.save()


def annotation_run_retry(annotation_run: AnnotationRun, upload_only=False) -> AnnotationRun:
    if upload_only and annotation_run.vcf_annotated_filename is None:
        msg = "Retry annotation run upload only requires annotation VCF to be written"
        raise ValueError(msg)

    annotation_range_lock = annotation_run.annotation_range_lock
    if annotation_range_lock is None:
        msg = "Can't retry annotation run with no annotation lock!"
        raise ValueError(msg)

    annotation_run.error_exception = None  # Clear so status won't be error
    annotation_run.task_id = None  # Allow celery jobs to get lock on it
    annotation_run.save()

    if upload_only:
        tasks = [
            delete_annotation_run_uploaded_data.si(annotation_run.pk),
        ]
    else:
        # Delete old AnnotationRun then try again.
        old_annotation_run = annotation_run
        annotation_run = AnnotationRun.objects.create()
        tasks = [
            delete_annotation_run.si(old_annotation_run.pk),
            assign_range_lock_to_annotation_run.si(annotation_run.pk, annotation_range_lock.pk),
        ]

    tasks.append(annotate_variants.si(annotation_run.pk))
    task = chain(tasks)
    task.apply_async()

    return annotation_run