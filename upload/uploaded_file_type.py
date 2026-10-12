import logging
from collections import defaultdict
from typing import Optional

from library.utils.file_utils import get_extension_without_gzip
from upload.import_task_factories.import_task_factory import (
    get_import_task_factories,
    get_import_task_factory_from_extension,
)
from upload.models import UploadData, UploadedVCF
from upload.tasks.vcf.genotype_vcf_tasks import reload_vcf_task
from upload.upload_processing import process_upload_pipeline


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


def _get_data_classes_by_file_type() -> dict[str, list[type[UploadData]]]:
    """ factories list their classes most specific first """
    data_classes_by_file_type = {}
    for itf in get_import_task_factories():
        file_type = itf.get_uploaded_file_type()
        data_classes_by_file_type[file_type] = itf.get_data_classes()
    return data_classes_by_file_type


def get_upload_data_for_uploaded_file(file_upload) -> Optional[UploadData]:
    """ The UploadData created by processing file_upload """
    classes = _get_data_classes_by_file_type().get(file_upload.file_type)
    if classes:
        for klazz in classes:
            try:
                return klazz.objects.get(file_upload=file_upload)
            except Exception:
                pass

    return None


def get_upload_data_by_file_upload_id(file_type_by_file_upload_id: dict[int, str]) -> dict[int, UploadData]:
    """ get_upload_data_for_uploaded_file for many files - a query per data class rather than per file """
    data_classes_by_file_type = _get_data_classes_by_file_type()
    classes_by_file_upload_id = {file_upload_id: data_classes_by_file_type.get(file_type, [])
                                 for file_upload_id, file_type in file_type_by_file_upload_id.items()}
    file_upload_ids_by_class = defaultdict(set)
    for file_upload_id, classes in classes_by_file_upload_id.items():
        for klazz in classes:
            file_upload_ids_by_class[klazz].add(file_upload_id)

    upload_data_by_class_and_id = {}
    for klazz, file_upload_ids in file_upload_ids_by_class.items():
        data_fields = [f.name for f in klazz._meta.concrete_fields if f.is_relation and f.name != "file_upload"]
        for upload_data in klazz.objects.filter(file_upload_id__in=file_upload_ids).select_related(*data_fields):
            upload_data_by_class_and_id.setdefault((klazz, upload_data.file_upload_id), upload_data)

    upload_data_by_file_upload_id = {}
    for file_upload_id, classes in classes_by_file_upload_id.items():
        for klazz in classes:
            if upload_data := upload_data_by_class_and_id.get((klazz, file_upload_id)):
                upload_data_by_file_upload_id[file_upload_id] = upload_data
                break
    return upload_data_by_file_upload_id


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
