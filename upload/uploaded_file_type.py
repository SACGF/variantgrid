import logging
from collections import defaultdict
from operator import itemgetter
from typing import Optional

from library.utils.file_utils import get_extension_without_gzip
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadData, UploadedVCF
from upload.tasks.vcf.genotype_vcf_tasks import reload_vcf_task
from upload.upload_processing import process_upload_pipeline


def get_import_tasks_by_extension():
    possible_tasks = defaultdict(list)
    for itf in get_import_task_factories():
        for ext in itf.get_possible_extensions():
            possible_tasks[ext].append(itf)
    return possible_tasks


def get_import_task_factory_from_extension(user, filename, file_extension):
    possible_tasks = get_import_tasks_by_extension()
    possible_for_extension = possible_tasks[file_extension]

    tasks = []
    for possible in possible_for_extension:
        processing_ability = possible.get_processing_ability(user, filename, file_extension)
        if processing_ability:
            tasks.append((int(processing_ability), possible))

    if tasks:
        logging.debug("tasks: %s", tasks)
        tasks = sorted(tasks, key=itemgetter(0), reverse=True)
        last_pa = None
        for pa, _ in tasks:
            if last_pa is not None:
                if pa == last_pa:
                    logging.warning("Task for extension %s had 2 processors with equal ability - can't decide!", file_extension)
                    return None
            else:
                last_pa = pa

        return tasks[0][1]

    logging.warning("No tasks found for %s", tasks)
    return None


def get_uploaded_file_type(file_upload, original_filename):
    """ When Django's UploadedFile saves, it may add extra chars to make it unique, eg:
        combined.vcf.gz => combined.vcf_HvzMe7j.gz
        So need to pass in original file name, which we'll use to get extension
    """
    filename = file_upload.get_filename()
    file_extension = get_extension_without_gzip(original_filename)
    import_task_factory = get_import_task_factory_from_extension(file_upload.user, filename, file_extension)
    if import_task_factory:
        return import_task_factory.get_uploaded_file_type()
    return None


def get_url_and_data_for_uploaded_file_data(file_upload):
    url = None
    upload_data = get_upload_data_for_uploaded_file(file_upload)
    if upload_data:
        url = upload_data.get_data_url()
    return url, upload_data


def get_upload_data_for_uploaded_file(file_upload) -> Optional[UploadData]:
    """ The UploadData created by processing file_upload - factories list their classes most specific first """
    data_classes_by_file_type = {}
    for itf in get_import_task_factories():
        file_type = itf.get_uploaded_file_type()
        data_classes_by_file_type[file_type] = itf.get_data_classes()

    classes = data_classes_by_file_type.get(file_upload.file_type)
    if classes:
        for klazz in classes:
            try:
                return klazz.objects.get(file_upload=file_upload)
            except Exception:
                pass

    return None


def reloads_vcf_in_place(upload_data) -> bool:
    """ True for every file type that loads a VCF - which is more than the '.vcf' ones, eg DRAGEN
        TSO500 AllFusions rows become a VCF (@see AbstractVCFImportTaskFactory).

        These reload through reload_vcf_task, which keeps the VCF and rebuilds its internal data.
        Deleting the UploadedVCF instead takes the VCF with it (@see pre_delete_uploaded_vcf), losing
        anything set on it by hand - eg a genome build the user picked because the file declared none. """
    return isinstance(upload_data, UploadedVCF)


def retry_upload_pipeline(upload_pipeline):
    upload_pipeline.remove_processing_files()

    logging.debug("retrying upload of %s", upload_pipeline)
    file_upload = upload_pipeline.file_upload

    upload_data = get_upload_data_for_uploaded_file(file_upload)
    if reloads_vcf_in_place(upload_data):
        task = reload_vcf_task.si(upload_pipeline.pk, upload_data.vcf_id)  # @UndefinedVariable
        task.apply_async()
    else:
        if upload_data and upload_data.created_by_pipeline:
            logging.debug("Type: %s, deleting file records: %s", file_upload.file_type, upload_data)
            upload_data.delete()

        # Re-use old UFPP so that it doesn't delete uploaded VCF
        upload_pipeline, *_ = process_upload_pipeline(upload_pipeline)
    return upload_pipeline
