import time
from datetime import datetime
from typing import Optional, Dict

import ijson
import requests
from dateutil import tz

from classification.models.evidence_key import EvidenceKeyMap
from classification.views.classification_view import BulkClassificationInserter
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector
from library.utils import make_json_safe_in_place, batch_iterator
from snpdb.models.models import Lab, Organization, Country
from sync.models.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.sync_runner import SyncRunner, register_sync_runner


@register_sync_runner(config={"type": {"shariant", "variantgrid"}, "direction": "download"})
class VariantGridDownloadSyncer(SyncRunner):

    def sync(self, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
        sync_destination = self.sync_destination

        config = sync_destination.config
        shariant = OAuthConnector.shariant_oauth_connector(sync_destination.sync_details)

        required_build = config.get('genome_build', 'GRCh37')
        params = {'share_level': 'public',
                  'type': 'json',
                  'build': required_build}

        exclude_labs = config.get('exclude_labs', None)
        if exclude_labs:
            params['exclude_labs'] = ','.join(exclude_labs)

        exclude_orgs = config.get('exclude_orgs', None)
        if exclude_orgs:
            params['exclude_orgs'] = ','.join(exclude_orgs)

        if not full_sync:
            last_download: SyncRun
            if last_download := SyncRun.objects.filter(destination=sync_destination, status=SyncStatus.SUCCESS).order_by('-created').first():
                if meta := last_download.meta:
                    if since_str := meta.get('server_date'):
                        since = datetime.strptime(since_str, "%a, %d %b %Y %H:%M:%S %Z").replace(tzinfo=tz.UTC).timestamp()
                        params['since'] = str(since)

        def sanitize_data(known_keys: EvidenceKeyMap, data: dict, source_url: str) -> dict:
            skipped_keys = []
            sanitized = {}
            for key, blob in data.items():
                if key == 'owner' or key == 'source_id':
                    pass
                elif known_keys.get(key).is_dummy:
                    skipped_keys.append(key)
                else:
                    sanitized[key] = blob

            source_url_note = None
            if skipped_keys:
                key_list = ', '.join(skipped_keys)
                source_url_note = f'Unable to import {key_list}'

            source_url_blob = {
                "value": source_url,
                "note": source_url_note
            }
            sanitized["source_url"] = source_url_blob

            return sanitized

        def shariant_download_to_upload(known_keys: EvidenceKeyMap, record: dict) -> Optional[Dict]:
            meta = record.get('meta', {})
            record_id = meta.get('id')

            source_url = shariant.url(f'classification/classification/{record_id}')
            data = record.get('data') or {}  # data=None for deleted records returned via 'since' which is patch-like
            data = sanitize_data(known_keys, data, source_url)
            data = {
                "id": record.get('id'),
                "publish": record.get('publish'),
                "data": data,
            }
            if delete := record.get('delete'):
                data['delete'] = delete
                return data

            lab_group_name = meta.get('lab_id')
            lab = Lab.objects.filter(group_name=lab_group_name).first()

            # in case we screwed up exclude, don't want to accidentally import over our own records with shariant copy
            if lab and not lab.external:
                return None

            if not lab:
                parts = lab_group_name.split('/')
                org, _ = Organization.objects.get_or_create(group_name=parts[0], defaults={"name": parts[0]})
                australia, _ = Country.objects.get_or_create(name='Australia')
                Lab.objects.create(
                    group_name=lab_group_name,
                    name=meta.get('lab_name'),
                    organization=org,
                    city='Unknown',
                    country=australia,
                    external=True,
                )
            return data

        run = SyncRun(destination=sync_destination, status=SyncStatus.IN_PROGRESS)
        run.save()
        try:
            response = requests.get(shariant.url('classification/api/classifications/export'),
                                    auth=shariant.auth(),
                                    params=params,
                                    stream=False,
                                    timeout=MINUTE_SECS)

            last_modified = response.headers.get('Last-Modified')
            evidence_keys = EvidenceKeyMap.instance()

            count = 0
            skipped = 0

            def records():
                nonlocal skipped
                for record in ijson.items(response.content, 'records.item'):
                    make_json_safe_in_place(record)
                    mapped = shariant_download_to_upload(evidence_keys, record)
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

            run.status = SyncStatus.SUCCESS
            run.meta = {
                "params": params,
                "records_upserted": count,
                "records_skipped": skipped,
                "server_date": last_modified
            }
            run.save()
        finally:
            if run.status == SyncStatus.IN_PROGRESS:
                run.status = SyncStatus.FAILED
                run.save()

        return run
