from functools import cached_property
from typing import Optional

import requests
from oauthlib.oauth2 import LegacyApplicationClient
from requests.auth import HTTPBasicAuth, AuthBase
from requests_oauthlib import OAuth2Session, OAuth2
from rest_framework.authentication import BaseAuthentication

from classification.enums import ShareLevel
from classification.views.exports import ClassificationExportFormatterMVL
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_mvl import FormatDetailsMVL, \
    FormatDetailsMVLFileFormat
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector
from snpdb.models import Lab, GenomeBuild, Organization
from sync.sync_runner import SyncRunner, register_sync_runner
import json

@register_sync_runner(config={"type": "alissa", "direction": "upload"})
class AlissaUploadSync(SyncRunner):

    def sync(self, full_sync: bool = False, max_rows: Optional[int] = None):

        exclude_group_names = self.get_config("exclude_sources")
        excludes = set()
        for source in exclude_group_names:
            if "/" in source:
                excludes.add(Lab.objects.filter(group_name=source).get())
            else:
                excludes.add(Organization.objects.filter(group_name=source).get())

        genome_build_name = self.get_config("genome_build")
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        format_details = FormatDetailsMVL()
        format_details.format = FormatDetailsMVLFileFormat.JSON

        mvl_id = int(self.get_config("mvl_id"))

        # only support full syncs for now
        exporter = ClassificationExportFormatterMVL(
            classification_filter=ClassificationFilter(
                user=admin_bot(),
                genome_build=genome_build,
                exclude_sources=excludes,
                min_share_level=ShareLevel.ALL_USERS
            ),
            format_details=format_details
        )
        alissa = OAuthConnector.shariant_oauth_connector(self.sync_destination.sync_details)
        auth = alissa.auth()
        for file in exporter.serve_in_memory():
            response = requests.post(
                url=alissa.url(f'managedvariantlists/{mvl_id}/import?curated=true&importOption=CONTRIBUTE'),
                auth=auth,
                json=json.loads(file),  # we're turning json into string to turn it back into json, probably a way we can send the already stringified version
                timeout=MINUTE_SECS,
            )
            response.raise_for_status()
