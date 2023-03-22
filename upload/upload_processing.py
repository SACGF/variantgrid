import logging
from typing import Tuple

from celery.result import AsyncResult

from eventlog.models import create_event
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadPipeline, ProcessingStatus, UploadedFileTypes, UploadedFile, UploadStepOrigin


def get_upload_processing_task(file_type, upload_pipeline):
    """ Finds subclass of ImportTaskFactory that returns file_type
        from get_uploaded_file_type() then uses it to create pipeline """
    factories = {}
    for itf in get_import_task_factories():
        uft = itf.get_uploaded_file_type()
        factories[uft] = itf

    factory = factories.get(file_type)
    if factory:
        import_task = factory.create_import_task(upload_pipeline)
        return import_task

    raise ValueError(f"No ImportTaskFactory found for file type '{file_type}'")


def process_uploaded_file(uploaded_file, run_async=True) -> Tuple[UploadPipeline, AsyncResult]:
    """ returns (UploadPipeline, result) """
    logging.debug("process_uploaded_file: filetype = %s", uploaded_file.file_type)

    try:
        upload_pipeline = UploadPipeline.objects.get(uploaded_file=uploaded_file)
        msg = f"There is already an UploadPipeline (pk={upload_pipeline.pk}) for uploaded file (pk='{uploaded_file.pk}')"
        raise ValueError(msg)
    except UploadPipeline.DoesNotExist:
        pass

    upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
    return process_upload_pipeline(upload_pipeline, run_async=run_async)


def process_upload_pipeline(upload_pipeline: UploadPipeline,
                            run_async=True) -> Tuple[UploadPipeline, AsyncResult]:
    """ Reuses the same upload_pipeline - relies on uploaded_file being set """

    if not upload_pipeline.uploaded_file:
        msg = "UploadPipeline.uploaded_file not set"
        raise ValueError(msg)

    # Only delete the steps that are in factory (will be recreated)
    # If anyone made any extra custom ones, keep them
    upload_steps = upload_pipeline.uploadstep_set
    upload_steps.filter(origin=UploadStepOrigin.IMPORT_TASK_FACTORY).delete()
    # Reset any other steps
    upload_steps.update(start_date=None, output_text=None, error_message='',
                        status=ProcessingStatus.CREATED,
                        items_to_process=0, items_processed=0, celery_task=None)
    upload_pipeline.status = ProcessingStatus.CREATED
    upload_pipeline.progress_status = "Created"
    upload_pipeline.items_to_process = 0
    upload_pipeline.items_processed = 0
    upload_pipeline.progress_percent = 0
    upload_pipeline.save()

    uploaded_file = upload_pipeline.uploaded_file
    user = uploaded_file.user
    create_event(user, f"upload_file_{uploaded_file.file_type}")

    task = get_upload_processing_task(uploaded_file.file_type, upload_pipeline)
    if task:
        if run_async:
            result = task.apply_async()
        else:
            result = task.apply()

        # Do as an update as other jobs may be modifying object
        UploadPipeline.objects.filter(id=upload_pipeline.pk).update(celery_task=result.id)
    else:
        result = None
    return upload_pipeline, result


def process_vcf_file(vcf_filename, name, user, import_source, run_async=True, file_type=UploadedFileTypes.VCF) -> Tuple[UploadPipeline, AsyncResult]:
    logging.info("process_vcf_file, path=%s", vcf_filename)
    uploaded_file = UploadedFile.objects.create(path=vcf_filename,
                                                import_source=import_source,
                                                name=name,
                                                user=user,
                                                file_type=file_type)
    ufpj, result = process_uploaded_file(uploaded_file, run_async)

    if not run_async:  # Should be done by now...
        ufpj = UploadPipeline.objects.get(pk=ufpj.pk)
    return ufpj, result

