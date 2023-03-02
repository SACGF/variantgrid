from functools import cached_property
from typing import Optional

import requests
from django.core.files.uploadedfile import InMemoryUploadedFile, SimpleUploadedFile
from oauthlib.oauth2 import LegacyApplicationClient
from requests.auth import AuthBase, HTTPBasicAuth
from requests_oauthlib import OAuth2Session, OAuth2

from classification.models import resolve_uploaded_url_to_handle, FileHandle, UploadedClassificationsUnmappedStatus, \
    UploadedClassificationsUnmapped
from classification.tasks.classification_import_map_and_insert_task import ClassificationImportMapInsertTask
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector
from snpdb.models import Lab
from sync.sync_runner import SyncRunner, register_sync_runner


@register_sync_runner(config={"type": "alissa", "direction": "download"})
class AlissaDownloadSync(SyncRunner):

    def sync(self, full_sync: bool = False, max_rows: Optional[int] = None):

        destination_lab = Lab.objects.filter(group_name=self.get_config("lab")).get()

        if not destination_lab.upload_location or not destination_lab.upload_automatic:
            raise ValueError("Can only download to labs with upload_location and marked as upload_automatic")

        mvl_id = int(self.get_config("mvl_id"))

        alissa = OAuthConnector.shariant_oauth_connector(self.sync_destination.sync_details)
        auth = alissa.auth()

        response = requests.get(
            url=alissa.url(f'managedvariantlists/{mvl_id}/export'),
            auth=auth,
            # we're turning json into string to turn it back into json, probably a way we can send the already stringified version
            timeout=MINUTE_SECS,
        )
        response.raise_for_status()

        if server_address := resolve_uploaded_url_to_handle(destination_lab.upload_location):
            sub_file = server_address.save(SimpleUploadedFile(self.sync_destination.name + ".json", content=response.content))

            status = UploadedClassificationsUnmappedStatus.Pending
            uploaded_file = UploadedClassificationsUnmapped.objects.create(
                url=sub_file.clean_url,
                filename=sub_file.filename,
                user=admin_bot(),
                lab=destination_lab,
                status=status,
                effective_modified=sub_file.modified,
                file_size=sub_file.size
            )

            task = ClassificationImportMapInsertTask.si(uploaded_file.pk)
            task.apply_async()
