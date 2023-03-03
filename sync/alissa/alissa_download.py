from django.core.files.uploadedfile import SimpleUploadedFile
from classification.models import resolve_uploaded_url_to_handle, UploadedClassificationsUnmappedStatus, \
    UploadedClassificationsUnmapped
from classification.tasks.classification_import_map_and_insert_task import ClassificationImportMapInsertTask
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from snpdb.models import Lab
from sync.sync_runner import SyncRunner, register_sync_runner, SyncRunInstance


@register_sync_runner(config={"type": "alissa", "direction": "download"})
class AlissaDownloadSyncer(SyncRunner):

    def sync(self, sync_run_instance: SyncRunInstance):
        if sync_run_instance.max_rows:
            raise ValueError("AlissaDownloadSyncer does not support max_rows")

        destination_lab = Lab.objects.filter(group_name=sync_run_instance.get_config("lab")).get()

        if not destination_lab.upload_location or not destination_lab.upload_automatic:
            raise ValueError("Can only download to labs with upload_location and marked as upload_automatic")

        mvl_id = int(sync_run_instance.get_config("mvl_id"))

        alissa = sync_run_instance.server_auth()

        sync_run_instance.run_start()

        response = alissa.get(
            url_suffix=f'managedvariantlists/{mvl_id}/export',
            timeout=MINUTE_SECS * 10,
        )
        response.raise_for_status()

        if server_address := resolve_uploaded_url_to_handle(destination_lab.upload_location):
            sub_file = server_address.save(SimpleUploadedFile(sync_run_instance.sync_destination.name + ".json", content=response.content))

            uploaded_file = UploadedClassificationsUnmapped.objects.create(
                url=sub_file.clean_url,
                filename=sub_file.filename,
                user=admin_bot(),
                lab=destination_lab,
                status=UploadedClassificationsUnmappedStatus.Pending,
                effective_modified=sub_file.modified,
                file_size=sub_file.size
            )

            task = ClassificationImportMapInsertTask.si(uploaded_file.pk)
            task.apply_async()

            sync_run_instance.run_completed(had_records=True)
