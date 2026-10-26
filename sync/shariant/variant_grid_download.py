import logging
import re
import time
from typing import Optional

import ijson

from classification.models.classification_inserter import BulkClassificationInserter
from classification.models.classification_ref import ClassificationRef
from classification.models.evidence_key import EvidenceKeyMap
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.log_utils import report_message
from library.oauth import ServerAuth
from library.utils import batch_iterator, make_json_safe_in_place
from snpdb.models.models import Country, Lab, Organization
from sync.models.models import SyncDestination
from sync.shariant.remote_response import request_json_with_retries
from sync.sync_runner import SyncRunInstance, SyncRunner, register_sync_runner

ORG_GROUP_NAME_PATTERN = re.compile(r'[\w\-]+')
LAB_GROUP_NAME_PATTERN = re.compile(r'[\w\-]+/[\w\-]+')


def _group_names_param(config: dict, key: str, pattern: re.Pattern) -> Optional[str]:
    """ Comma separated group names for the export API, which silently ignores a name it can't match -
        so fail on a malformed one rather than download records we meant to exclude """
    if group_names := config.get(key):
        for group_name in group_names:
            if not (isinstance(group_name, str) and pattern.fullmatch(group_name)):
                raise ValueError(f"SyncDestination config {key} has malformed group name {group_name!r}")
        return ','.join(group_names)
    return None


def _sanitize_data(known_keys: EvidenceKeyMap, data: dict, source_url: str) -> dict:
    skipped_keys = []
    sanitized = {}
    for key, blob in data.items():
        if key in ('owner', 'source_id'):
            pass
        elif known_keys.get(key).is_dummy:
            skipped_keys.append(key)
        else:
            sanitized[key] = blob

    source_url_note = None
    if skipped_keys:
        key_list = ', '.join(skipped_keys)
        source_url_note = f'Unable to import {key_list}'

    sanitized["source_url"] = {
        "value": source_url,
        "note": source_url_note
    }
    return sanitized


def _record_lab_group_name(record_id) -> Optional[str]:
    """ The lab the inserter files a record under, from its "org/lab/lab_record_id" id - None if the id isn't one """
    if isinstance(record_id, str):
        lab_group_name, lab_record_id, _version = ClassificationRef.parse_id_str(record_id)
        if lab_group_name and lab_record_id and LAB_GROUP_NAME_PATTERN.fullmatch(lab_group_name):
            return lab_group_name
    return None


def _create_external_lab(sync_destination: SyncDestination, lab_group_name: str, lab_name) -> Lab:
    org_group_name, lab_short_group_name = lab_group_name.split('/')
    if not (isinstance(lab_name, str) and lab_name.strip()):
        lab_name = lab_short_group_name
    org, org_created = Organization.objects.get_or_create(group_name=org_group_name, defaults={"name": org_group_name})
    australia, _ = Country.objects.get_or_create(name='Australia')
    lab = Lab.objects.create(
        group_name=lab_group_name,
        name=lab_name,
        organization=org,
        city='Unknown',
        country=australia,
        external=True,
    )
    report_message("Sync download created external lab", extra_data={
        "target": lab_group_name,
        "lab_name": lab_name,
        "organization_created": org_created,
        "sync_destination": sync_destination.name,
    })
    return lab


def shariant_download_to_upload(server_auth: ServerAuth, sync_destination: SyncDestination,
                                known_keys: EvidenceKeyMap, record_json: dict) -> Optional[dict]:
    """ A remote export record as a BulkClassificationInserter record, or None to skip it.
        The inserter (as the admin bot) files the record under the lab named in its id, so that lab - for deletes
        too - must be external or not yet exist: in case exclude_labs is wrong, never import over our own records """
    record_id = record_json.get('id')
    lab_group_name = _record_lab_group_name(record_id)
    if not lab_group_name:
        # the remote names the Lab and Organization we create, so hold it to the group name format
        report_message("Sync download skipped record with malformed id", extra_data={
            "target": repr(record_id),
            "sync_destination": sync_destination.name,
        })
        return None

    lab = Lab.objects.filter(group_name=lab_group_name).first()
    if lab and not lab.external:
        return None

    meta = record_json.get('meta') or {}
    delete = record_json.get('delete')
    if not lab:
        if delete:
            return None  # no records of an unknown lab are held here to withdraw
        _create_external_lab(sync_destination, lab_group_name, meta.get('lab_name'))

    source_url = server_auth.url(f"classification/classification/{meta.get('id')}")
    data = record_json.get('data') or {}  # data=None for deleted records returned via 'since' which is patch-like
    upload = {
        "id": record_id,
        "publish": record_json.get('publish'),
        "data": _sanitize_data(known_keys, data, source_url),
    }
    if delete:
        upload['delete'] = delete
    return upload


@register_sync_runner(config={"type": {"shariant", "variantgrid"}, "direction": "download"})
class VariantGridDownloadSyncer(SyncRunner):

    def sync(self, sync_run_instance: SyncRunInstance):
        if sync_run_instance.max_rows:
            raise ValueError("VariantGridDownloadSyncer does not support max_rows")

        sync_destination = sync_run_instance.sync_destination

        config = sync_destination.config
        other_variant_grid = ServerAuth.for_sync_details(sync_destination.sync_details)

        required_build = config.get('genome_build', 'GRCh37')
        params = {'share_level': 'public',
                  'type': 'json',
                  'build': required_build}

        if exclude_labs := _group_names_param(config, 'exclude_labs', LAB_GROUP_NAME_PATTERN):
            params['exclude_labs'] = exclude_labs

        if exclude_orgs := _group_names_param(config, 'exclude_orgs', ORG_GROUP_NAME_PATTERN):
            params['exclude_orgs'] = exclude_orgs

        if not sync_run_instance.full_sync:
            if since := sync_run_instance.last_success_server_date():
                params['since'] = str(since.timestamp())

        sync_run_instance.run_start()
        download_start = time.time()
        # stream=False so a body cut short in transit fails (and retries) here, before any record is inserted
        response = request_json_with_retries(
            lambda: other_variant_grid.get(
                url_suffix='classification/api/classifications/export',
                params=params,
                stream=False,
                timeout=MINUTE_SECS
            ),
            description=f"{sync_destination}: download export",
        )
        logging.info("%s: downloaded export in %.1fs", sync_destination, time.time() - download_start)

        last_modified = response.headers.get('Last-Modified')
        evidence_keys = EvidenceKeyMap.instance()

        count = 0
        skipped = 0

        def records():
            nonlocal skipped
            for record_item in ijson.items(response.content, 'records.item'):
                make_json_safe_in_place(record_item)
                mapped = shariant_download_to_upload(other_variant_grid, sync_destination, evidence_keys, record_item)
                if mapped:
                    yield mapped
                else:
                    skipped = skipped + 1

        first_batch = True
        for batch in batch_iterator(records(), batch_size=50):
            inserter = BulkClassificationInserter(user=admin_bot(), force_publish=True)
            try:
                # give the process 10 seconds to breath between batches of 50 classifications
                # in case we're downloading giant chunks of data
                if not first_batch:
                    time.sleep(10)
                first_batch = False

                for record in batch:
                    inserter.insert(record)
                    count = count + 1
            finally:
                inserter.finish()
            logging.info("%s: upserted %d record(s) so far (%d skipped)", sync_destination, count, skipped)

        sync_run_instance.run_completed(
            had_records=bool(count),
            meta={
                "params": params,
                "records_upserted": count,
                "records_skipped": skipped,
                "server_date": last_modified
            }
        )
