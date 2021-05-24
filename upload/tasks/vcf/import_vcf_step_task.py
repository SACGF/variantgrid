from celery.app.task import Task
from celery.canvas import chord, chain
from django.db.models.aggregates import Min, Max
from django.db.models.expressions import F
from django.utils import timezone
import celery
import logging
import subprocess

from library.log_utils import get_traceback
from library.utils import import_class
from upload.models import UploadPipeline, ProcessingStatus, \
    UploadStep, PipelineFailedJobTerminateEarlyException, VCFPipelineStage, SkipUploadStepException


class ImportVCFStepTask(Task):

    def process_items(self, upload_step):
        raise NotImplementedError()

    @staticmethod
    def check_pipeline_stage(upload_pipeline, upload_step):
        """ If you're the last step of your pipeline_stage, kick off any dependencies """

        if upload_step.pipeline_stage_dependency == VCFPipelineStage.FINISH:
            return  # Nothing to kick off

        other_unfinished_upload_steps_qs = upload_pipeline.uploadstep_set.filter(end_date__isnull=True).exclude(pk=upload_step.pk)

        pipeline_stage = upload_step.pipeline_stage
        if pipeline_stage:
            other_unfinished_upload_steps_for_stage_qs = other_unfinished_upload_steps_qs.filter(pipeline_stage=pipeline_stage)
            waiting_upload_steps_qs = other_unfinished_upload_steps_qs.filter(pipeline_stage_dependency=pipeline_stage,
                                                                              status=ProcessingStatus.CREATED,
                                                                              start_date__isnull=True)

            if not other_unfinished_upload_steps_for_stage_qs.exists() and waiting_upload_steps_qs.exists():
                schedule_pipeline_stage_steps.apply_async((upload_pipeline.pk, pipeline_stage))
                return  # don't check finish jobs below

        # Finish jobs - when all steps other than those depending on FINISH are done
        non_finish_jobs_qs = other_unfinished_upload_steps_qs.exclude(pipeline_stage_dependency=VCFPipelineStage.FINISH)
        if not non_finish_jobs_qs.exists():
            schedule_pipeline_stage_steps.apply_async((upload_pipeline.pk, VCFPipelineStage.FINISH))

    def run(self, upload_step_id, pipeline_fraction=0):

        def load_upload_step():
            return UploadStep.objects.get(pk=upload_step_id)

        upload_step = load_upload_step()
        upload_pipeline = upload_step.upload_pipeline
        upload_pipeline_qs = UploadPipeline.objects.filter(id=upload_pipeline.pk)

        try:
            this_task = self.request.id
            if upload_step.celery_task is not None and upload_step.celery_task != this_task:
                params = (upload_step.upload_pipeline.pk, upload_step.pk, this_task, upload_step.celery_task)
                message = "It appears Pipeline %d, step %d is being run twice - this task='%s', upload_step.celery_task already set to '%s'" % params
                raise Exception(message)

            upload_step.start()
            upload_step.celery_task = this_task
            upload_step.save()

            if upload_pipeline.status != ProcessingStatus.PROCESSING:
                upload_step.status = ProcessingStatus.SKIPPED
            else:
                items_processed = self.process_items(upload_step)
                upload_step = load_upload_step()  # Reload from DB (other processes may have modified)
                upload_step.items_processed = items_processed

                # If expected and actual items processed are provided, check they match
                if upload_step.items_to_process is not None and items_processed is not None:
                    if upload_step.items_to_process != upload_step.items_processed:
                        params = (upload_step.pk, upload_step.items_to_process, upload_step.items_processed)
                        msg = "UploadStep (%d) has items_to_process (%d) != items_processed (%d)" % params
                        raise ValueError(msg)

                upload_step.status = ProcessingStatus.SUCCESS
                if pipeline_fraction:
                    upload_pipeline_qs.update(progress_percent=F('progress_percent') + pipeline_fraction)

        except PipelineFailedJobTerminateEarlyException:
            logging.error("ImportVCFStepTask: Got PipelineFailedJobTerminateEarlyException")
            upload_step.status = ProcessingStatus.TERMINATED_EARLY
        except SkipUploadStepException:
            upload_step.status = ProcessingStatus.SKIPPED
        except subprocess.CalledProcessError as e:
            message = e.output
            upload_step.error_exception(message)
            upload_pipeline.error(message)
        except Exception as e:
            message = str(e)
            if not upload_pipeline_qs.exists():
                logging.warning("UploadPipeline was deleted, causing error:")
                logging.warning(message)
                return

            upload_step.error_exception(message)
            upload_pipeline.error(message)

        upload_step.end_date = timezone.now()
        upload_step.save()

        upload_step.close_sub_steps()

        # Do this after this step is closed, so that if there is a race condition - it'll be to trigger
        # schedule_pipeline_stage_steps >1 times instead of 0 (will be handled ok as it's in a single worker queue)
        upload_pipeline = upload_pipeline_qs.get()  # Reload from DB
        if upload_pipeline.status == ProcessingStatus.PROCESSING:
            self.check_pipeline_stage(upload_pipeline, upload_step)


@celery.task
def schedule_pipeline_stage_steps(upload_pipeline_id, pipeline_stage):
    """ Only want to do this once per VCF import, so run in own task on single queue to avoid race conditions.
        May be executed multiple times but will only run jobs the 1st time """

    upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)

    try:
        logging.info("schedule_pipeline_stage_steps: %s", pipeline_stage)
        parallel_tasks = []
        pipeline_stage_dependent_steps_qs = upload_pipeline.uploadstep_set.filter(pipeline_stage_dependency=pipeline_stage,
                                                                                  start_date__isnull=True)
        for upload_step in pipeline_stage_dependent_steps_qs.order_by("sort_order"):
            # Set start date so that they will not be loaded in pipeline_stage_dependent_steps_qs, were it to run again.
            upload_step.start_date = timezone.now()
            upload_step.save()

            task_class = import_class(upload_step.script)
            parallel_tasks.append(task_class.si(upload_step.pk, 0))

        if pipeline_stage == VCFPipelineStage.FINISH:
            # Getting infinite celery.chord_unlock - see #175
            # chord_task = chord(parallel_tasks, pipeline_success_task.si(upload_pipeline_id))
            # So using chain to run in sequence not parallel

            # Add on pipeline success task at the end of FINISH jobs
            parallel_tasks.append(pipeline_success_task.si(upload_pipeline_id))
            chain_task = chain(parallel_tasks)
            chain_task.apply_async()
        else:
            for task in parallel_tasks:
                task.apply_async()
    except AttributeError:
        message = get_traceback()
        message += "\nYou probably didn't register your Celery class via app.register_task().\n"
        upload_pipeline.error(message)
    except:
        message = get_traceback()
        upload_pipeline.error(message)


@celery.task
def pipeline_start_task(upload_pipeline_id):
    upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)
    upload_pipeline.start()


@celery.task
def pipeline_success_task(upload_pipeline_id):
    upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)
    if upload_pipeline.status == ProcessingStatus.PROCESSING:
        steps = upload_pipeline.uploadstep_set.all()
        values = steps.aggregate(pipeline_start=Min('start_date'), pipeline_end=Max('end_date'))
        pipeline_start = values['pipeline_start']
        pipeline_end = values['pipeline_end']
        processing_seconds_wall_time = (pipeline_end - pipeline_start).total_seconds()
        try:
            processing_seconds_cpu_time = sum((s.end_date - s.start_date).total_seconds() for s in steps)
        except TypeError:
            processing_seconds_cpu_time = None

        upload_pipeline.success(processing_seconds_wall_time=processing_seconds_wall_time,
                                processing_seconds_cpu_time=processing_seconds_cpu_time)
