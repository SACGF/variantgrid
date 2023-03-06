from datetime import datetime
from enum import Enum
from typing import Optional

from classification.enums import ShareLevel
from classification.views.exports import ClassificationExportFormatterMVL
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_mvl import FormatDetailsMVL, \
    FormatDetailsMVLFileFormat
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from snpdb.models import Lab, GenomeBuild, Organization
from sync.sync_runner import SyncRunner, register_sync_runner, SyncRunInstance
import json


class AlissaImportOption(str, Enum):
    CONTRIBUTE = "CONTRIBUTE"
    MIRROR = "MIRROR"


@register_sync_runner(config={"type": "alissa", "direction": "upload"})
class AlissaUploadSyncer(SyncRunner):
    """
    Alissa uploads
    required config parameters:
        mvl_id: Index of the mvl to be uploaded to from this VariantGrid instance
        genome_build: The genome build to export using
        exclude_sources: An array of group names (either org or lab) to exclude from upload,
            should include the lab that we're uploading to.
    """

    def sync(self, sync_run_instance: SyncRunInstance):
        if sync_run_instance.max_rows:
            raise ValueError("AlissaUploadSyncer does not support max_rows")

        exclude_group_names = sync_run_instance.get_config("exclude_sources")
        excludes = set()
        for source in exclude_group_names:
            if "/" in source:
                excludes.add(Lab.objects.filter(group_name=source).get())
            else:
                excludes.add(Organization.objects.filter(group_name=source).get())

        genome_build_name = sync_run_instance.get_config("genome_build")
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        format_details = FormatDetailsMVL()
        format_details.format = FormatDetailsMVLFileFormat.JSON

        mvl_id = int(sync_run_instance.get_config("mvl_id"))

        since: Optional[datetime] = None
        if not sync_run_instance.full_sync:
            since = sync_run_instance.last_success_server_date()

        exporter = ClassificationExportFormatterMVL(
            classification_filter=ClassificationFilter(
                user=admin_bot(),
                genome_build=genome_build,
                exclude_sources=excludes,
                min_share_level=ShareLevel.ALL_USERS,
                since=since
            ),
            format_details=format_details
        )

        alissa = sync_run_instance.server_auth()
        uploaded_any_rows = False

        # For full syncs, set the ImportOption to "MIRROR", which will delete everything and replace it with our first post
        # then subsequent posts will just add
        import_option = AlissaImportOption.MIRROR if sync_run_instance.full_sync else AlissaImportOption.CONTRIBUTE

        for file in exporter.serve_in_memory():
            if exporter.row_count > 0:
                uploaded_any_rows = True
                sync_run_instance.run_start()
                params = {
                    "curated": "true",
                    "importOption": import_option.value
                }
                import_option = AlissaImportOption.CONTRIBUTE

                response = alissa.post(
                    url_suffix=f'managedvariantlists/{mvl_id}/import',
                    params=params,
                    json=json.loads(file),  # we're turning json into string to turn it back into json, probably a way we can send the already stringified version
                    timeout=MINUTE_SECS,
                )
                response.raise_for_status()

        sync_run_instance.run_completed(
            had_records=uploaded_any_rows,
            meta={
                "server_date": exporter.classification_filter.last_modified_header,
                "rows_uploaded": exporter.row_count
            }
        )
