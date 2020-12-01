from collections import defaultdict
import logging
from operator import itemgetter

from library.file_utils import get_extension_without_gzip
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadedFileTypes
from upload.tasks.vcf.genotype_vcf_tasks import reload_vcf_task
from upload.upload_processing import process_upload_pipeline


def get_import_task_factory_from_extension(user, filename, file_extension):
    possible_tasks = defaultdict(list)
    for itf in get_import_task_factories():
        for ext in itf.get_possible_extensions():
            possible_tasks[ext].append(itf)

    possible_for_extension = possible_tasks[file_extension]

    tasks = []
    for possible in possible_for_extension:
        processing_ability = possible.get_processing_ability(user, filename, file_extension)
        if processing_ability:
            tasks.append((int(processing_ability), possible))

    if tasks:
        print(f"tasks: {tasks}")
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


def get_uploaded_file_type(uploaded_file, original_filename):
    """ When Django UploadedFile saves, it may add extra chars to make it unique, eg:
        combined.vcf.gz => combined.vcf_HvzMe7j.gz
        So need to pass in original file name, which we'll use to get extension
    """
    filename = uploaded_file.get_filename()
    file_extension = get_extension_without_gzip(original_filename)
    import_task_factory = get_import_task_factory_from_extension(uploaded_file.user, filename, file_extension)
    if import_task_factory:
        return import_task_factory.get_uploaded_file_type()
    return None


def get_url_and_data_for_uploaded_file_data(uploaded_file):
    url = None
    upload_data = get_upload_data_for_uploaded_file(uploaded_file)
    if upload_data:
        try:
            data = upload_data.get_data()
            url = data.get_absolute_url()
        except:
            pass
    return url, upload_data


def get_upload_data_for_uploaded_file(uploaded_file):
    uploaded_file_classes = {}
    for itf in get_import_task_factories():
        file_type = itf.get_uploaded_file_type()
        uploaded_file_classes[file_type] = itf.get_data_classes()

    classes = uploaded_file_classes.get(uploaded_file.file_type)
    if classes:
        for klazz in classes:
            try:
                return klazz.objects.get(uploaded_file=uploaded_file)
            except Exception:
                pass

    return None


def retry_upload_pipeline(upload_pipeline):
    upload_pipeline.remove_processing_files()

    logging.debug("retrying upload of %s", upload_pipeline)
    uploaded_file = upload_pipeline.uploaded_file

    upload_data = get_upload_data_for_uploaded_file(uploaded_file)
    if upload_pipeline.file_type == UploadedFileTypes.VCF:
        if upload_data.vcf:
            vcf_id = upload_data.vcf.pk
        else:
            vcf_id = None
        task = reload_vcf_task.si(upload_pipeline.pk, vcf_id)  # @UndefinedVariable
        task.apply_async()
    else:
        if upload_data:
            logging.debug("Type: %s, deleting file records: %s", uploaded_file.file_type, upload_data)
            upload_data.delete()

        # Re-use old UFPP so that it doesn't delete uploaded VCF
        upload_pipeline, *_ = process_upload_pipeline(upload_pipeline)
    return upload_pipeline
