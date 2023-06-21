from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Optional, Set

from classification.enums import ShareLevel
from classification.views.exports import ClassificationExportFormatterMVL
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_mvl import FormatDetailsMVL, \
    FormatDetailsMVLFileFormat
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.log_utils import AdminNotificationBuilder, report_exc_info
from library.utils import ExportRow, export_column
from snpdb.models import Lab, GenomeBuild, Organization
from sync.models import SyncRun
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

        excludes = set()
        if exclude_group_names := sync_run_instance.get_config("exclude_sources", mandatory=False):
            for source in exclude_group_names:
                if "/" in source:
                    excludes.add(Lab.objects.filter(group_name=source).get())
                else:
                    excludes.add(Organization.objects.filter(group_name=source).get())

        # this will only be used for testing, note that it's mutually exclusive with exclude
        includes = set()
        if include_group_names := sync_run_instance.get_config("include_sources", mandatory=False):
            for source in include_group_names:
                includes.add(Lab.objects.filter(group_name=source).get())

        if not excludes and not includes:
            raise ValueError("Either exclude_sources or include_sources must be provided")

        genome_build_name = sync_run_instance.get_config("genome_build")
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        format_details = FormatDetailsMVL()
        format_details.format = FormatDetailsMVLFileFormat.JSON
        if mapping := sync_run_instance.get_config("classification_mapping"):
            classification_mapping = {
                'B': 'BENIGN',
                'LB': 'LIKELY_BENIGN',
                'VUS': 'VOUS',
                'LP': 'LIKELY_PATHOGENIC',
                'P': 'PATHOGENIC'
            }
            classification_mapping.update(mapping)
            format_details.classification_mapping = classification_mapping

        mvl_id = int(sync_run_instance.get_config("mvl_id"))

        since: Optional[datetime] = None
        if not sync_run_instance.full_sync:
            since = sync_run_instance.last_success_server_date()

        exporter = ClassificationExportFormatterMVL(
            classification_filter=ClassificationFilter(
                user=admin_bot(),
                genome_build=genome_build,
                exclude_sources=excludes,
                include_sources=includes,
                min_share_level=ShareLevel.ALL_USERS,
                since=since,
                rows_per_file=1000, # Alissa can handle 10,000 but 1,000 at a time just seems safer
                row_limit=sync_run_instance.max_rows
            ),
            format_details=format_details
        )

        alissa = sync_run_instance.server_auth()
        uploaded_any_rows = False

        # For full syncs, set the ImportOption to "MIRROR", which will delete everything and replace it with our first post
        # then subsequent posts will just add
        partial_upload = (not sync_run_instance.full_sync) or sync_run_instance.max_rows
        import_option = AlissaImportOption.CONTRIBUTE if partial_upload else AlissaImportOption.MIRROR
        response_jsons = []

        total_failed = 0
        total_differs = 0
        total_imported = 0
        for file in exporter.serve_in_memory():
            if exporter.row_count > 0:
                json_data = json.loads(file)

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
                    json=json_data,  # we're turning json into string to turn it back into json, probably a way we can send the already stringified version
                    timeout=MINUTE_SECS,
                )
                response.raise_for_status()

                try:
                    if response_json := response.json():
                        total_failed += int(response_json.get("numberFailed"))
                        total_differs += int(response_json.get("numberDiffers"))
                        total_imported += int(response_json.get("numberImported"))

                        response_jsons.append(response_json)
                        if response_error := response_json.get("error"):
                            notify = AdminNotificationBuilder(message="Error Uploading")
                            notify.add_field("Sync Destination", sync_run_instance.name)
                            notify.add_field("Error", response_error)
                            notify.send()
                        elif numberFailed := int(response_json.get("numberFailed")):
                            if numberFailed > 0:
                                notify = AdminNotificationBuilder(message="Error Uploading")
                                notify.add_field("Sync Destination", sync_run_instance.sync_destination.name)
                                notify.add_field("Failures", numberFailed)

                                failure: str
                                for failure in response_json.get("failures"):
                                    if "\t" in failure:
                                        parts = failure.split("\t")
                                        message = parts[0]
                                        json_data_str = parts[1]
                                        try:
                                            json_data = json.loads(json_data_str)
                                            transcript = json_data.get("transcript")
                                            c_nomen = json_data.get("cNomen")
                                            notify.add_markdown(f"{transcript}:{c_nomen} - \"{message}\"")
                                        except ValueError:
                                            notify.add_markdown(failure)
                                    else:
                                        notify.add_markdown(failure)

                                notify.send()
                except:
                    report_exc_info()
                    pass

        since_timestamp = None
        if since:
            since_timestamp = int(since.timestamp())
        sync_run_instance.run_completed(
            had_records=uploaded_any_rows,
            meta={
                "since": since_timestamp,
                "server_date": exporter.classification_filter.last_modified_header,
                "rows_sent": exporter.row_count,
                "total_failed": total_failed,
                "total_differs": total_differs,
                "total_imported": total_imported,
                "responses": response_jsons
            }
        )

    def report_on(self, sync_run: SyncRun):

        @dataclass(frozen=True)
        class AlissaRowInfoExport(ExportRow):
            c_hgvs: str
            severity: str
            message: str

            @export_column("c.HGVS")
            def _c_hgvs(self):
                return self.c_hgvs

            @export_column("severity")
            def _severity(self):
                return self.severity

            @export_column("message")
            def _message(self):
                return self.message

        all_issues: Set[AlissaRowInfoExport] = set()

        for response_dict in sync_run.meta.get("responses", {}):
            for info_line in response_dict.get("infos", []):
                info_parts = [p.strip().replace(" ", ":") for p in info_line.split("\t")]
                all_issues.add(AlissaRowInfoExport(c_hgvs=info_parts[0], severity="info", message=info_parts[1]))

            for failure in response_dict.get("failures"):
                parts = failure.split("\t")
                message = parts[0]
                json_data = json.loads(parts[1])

                transcript = json_data.get("transcript")
                c_nomen = json_data.get("cNomen")
                all_issues.add(
                    AlissaRowInfoExport(
                        c_hgvs=f"{transcript}:{c_nomen}",
                        severity="failure",
                        message=message
                    )
                )

        sorted_data = list(sorted(all_issues, key=lambda x: (x.severity, x.c_hgvs)))
        return AlissaRowInfoExport.streaming_csv(data=sorted_data, filename=f"sync_run_{sync_run.pk}")
