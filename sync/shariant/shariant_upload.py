import socket
from typing import Dict, Optional, Iterable, List, TypeVar, Union, Tuple

import requests
from django.db.models import QuerySet

from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector

from sync.models.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter
from classification.enums.classification_enums import ShareLevel
from classification.models import EvidenceKeyMap
from classification.models.classification import ClassificationModification
from classification.models.classification_utils import ClassificationJsonParams
from sync.shariant.query_json_filter import QueryJsonFilter

# add variant_type to private fields as the key has been deprecated
SHARIANT_PRIVATE_FIELDS = [
    'age_units', 'dob', 'family_id', 'internal_use', 'patient_id', 'patient_summary', 'sample_id', 'variant_type'
]


def insert_nones(data: Dict) -> Dict:
    for e_key in EvidenceKeyMap.instance().all_keys:
        key = e_key.key
        if key not in data:
            data[key] = None
    return data


T = TypeVar("T")


def batch_iterator_end(iterable: Iterable[T], batch_size: int = 10) -> Iterable[Union[List[T], bool]]:
    """
    Creates an iterator of list of T from an iterator of T, as well as providing a boolean to indicate if
    this is the final batch
    :param iterable: Iterable of T
    :param batch_size: Max number of items allowed in a batch (all but the last batch should be this size)
    :return: Union of list of T, and a boolean True if the batch is final, False otherwise
    """
    batch = []
    pending_batch = None
    for record in iterable:
        if pending_batch:
            yield pending_batch, False
            pending_batch = None

        batch.append(record)
        if len(batch) >= batch_size:
            pending_batch = batch
            batch = []

    if pending_batch:
        yield pending_batch, True
    if batch:
        yield batch, True


class ClassificationUploader:

    def __init__(self, sync_destination: SyncDestination):
        self.sync_destination = sync_destination

        config = self.sync_destination.config
        self.filters = config.get('filters', {})
        mapping = config.get('mapping', {})

        config_key = config["sync_details"]
        self.shariant = OAuthConnector.shariant_oauth_connector(config_key)
        self.lab_mappings = mapping.get('labs', {})
        self.share_level_mappings = mapping.get('share_levels', {})
        self.user_mappings = mapping.get('users', {})
        self.historical_converter = HistoricalEKeyConverter()

    def records_to_sync(self, apply_filters: bool = True, full_sync: bool = False) -> QuerySet[ClassificationModification]:
        qs = ClassificationModification.objects.filter(
            is_last_published=True,
            share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
            classification__lab__group_name__in=self.lab_mappings.keys()
        )
        if apply_filters and self.filters:
            q = QueryJsonFilter.classification_value_filter().convert_to_q(self.filters)
            print(q)
            qs = qs.filter(q)

        if not full_sync:
            qs = ClassificationModificationSyncRecord.filter_out_synced(
                qs=qs,
                destination=self.sync_destination
            )

        return qs

    def classification_to_json(self, vcm: ClassificationModification) -> Dict:
        raw_json = vcm.as_json(ClassificationJsonParams(current_user=admin_bot(), include_data=True, include_messages=False, strip_complicated=True, api_version=2))
        formatted_json = {}

        lab = vcm.classification.lab
        mapped_lab_name = self.lab_mappings.get(lab.group_name)
        if mapped_lab_name is True:
            mapped_lab_name = lab.group_name

        share_level = vcm.share_level
        share_level = self.share_level_mappings.get(share_level, share_level)

        # might need to map the lab group names if we don't map them all into shariant
        formatted_json['id'] = mapped_lab_name + '/' + vcm.classification.lab_record_id
        data = raw_json.get('data')
        for dont_share in SHARIANT_PRIVATE_FIELDS:
            data.pop(dont_share, None)

        # inject AMP if HORAK is used and vice versa (since VG3 only has a single value for assertion method)
        if data.get("assertion_method", {}).get("value") == "horak" or data.get("assertion_method", {}).get("value") == "amp":
            data["assertion_method"] = {"value": ["horak", "amp"]}

        # no need to screw around with owner
        user = data.get('owner', {}).get('value')
        if not user:
            user = vcm.user.username

        user = self.user_mappings.get(user, user)
        data['owner'] = {'value': user}

        data = insert_nones(data)
        data = self.historical_converter.to_shariant(vcm, data)
        formatted_json['overwrite'] = data
        formatted_json['publish'] = share_level
        # only allow for withdrawing, doesn't allow for un-withdrawing
        # but Shariant should get a notification that a withdrawn classification is being re-sent
        if raw_json.get('withdrawn', False):
            formatted_json['delete'] = True

        return formatted_json

    def sync(self, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
        qs = self.records_to_sync(full_sync=full_sync)

        rows_uploaded = 0
        run = SyncRun(destination=self.sync_destination, status=SyncStatus.IN_PROGRESS)
        run.save()

        if not qs.exists():
            run.status = SyncStatus.NO_RECORDS
            run.save()
        else:
            if max_rows is not None:
                qs = qs[:max_rows]

            site_name = socket.gethostname().lower().split('.')[0].replace('-', '')
            try:
                for batch, finished in batch_iterator_end(qs, batch_size=50):
                    # providing import_id and status:complete when we're done lets Shariant know
                    # when the upload has been completed
                    json_to_send = {
                        "records": [self.classification_to_json(vcm) for vcm in batch],
                        "import_id": site_name
                    }
                    if finished:
                        json_to_send["status"] = "complete"
                    # print(json.dumps(json_to_send))
                    auth = self.shariant.auth()

                    response = requests.post(
                        self.shariant.url('classification/api/classifications/v2/record/'),
                        auth=auth,
                        json=json_to_send
                    )
                    response.raise_for_status()
                    # results are sent back in an array in the same order they were sent up
                    results = response.json().get('results')

                    for record, result in zip(batch, results):
                        rows_uploaded += 1
                        ClassificationModificationSyncRecord.objects.create(run=run, classification_modification=record,
                                                                            meta=result)

                run.status = SyncStatus.SUCCESS
                run.save()
            finally:
                # a bit of a hack so I can let the exception just fall through to rollbar
                # while still updating the status of the run
                if run.status == SyncStatus.IN_PROGRESS:
                    run.status = SyncStatus.FAILED
                    run.save()

        return run

    def get_sync_match_status_and_remote_url(self, classification) -> Tuple[bool, Optional[bool], Optional[str]]:
        """
            :return:
                * sync_match - True/False if sync will collect classification or filter it out
                * sync_success - status of last sync run
                * remote_url - Url on remote system
        """
        qs = self.records_to_sync(full_sync=True).filter(classification=classification)
        sync_match = False
        sync_success = None
        remote_url = None
        if cm := qs.first():  # Should only be 1
            sync_match = True
            if cm_sync_record := ClassificationModificationSyncRecord.objects.filter(run__destination=self.sync_destination,
                                                                                    classification_modification=cm).latest():
                if sync_success := cm_sync_record.success:
                    # Search is a bit of a hack as we don't have a way to access it directly
                    remote_url = self.shariant.url(f'variantopedia/search?search={classification.lab_record_id}')

        return sync_match, sync_success, remote_url


# method is here for backwards compatibility, prefer to use ClassificationUploader directly
def sync_shariant_upload(sync_destination: SyncDestination, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
    return ClassificationUploader(sync_destination).sync(full_sync=full_sync, max_rows=max_rows)


def shariant_upload_status(classification):
    destination_status_and_url = {}
    sync_enabled = True
    # TODO: Remove this code after testing
    if socket.gethostname() == "sapathtest":
        sync_enabled = False  # Special case testing code

    for sd in SyncDestination.objects.filter(config__direction='upload', enabled=sync_enabled):
        uploader = ClassificationUploader(sd)
        sync_match, sync_success, remote_url = uploader.get_sync_match_status_and_remote_url(classification)
        if sync_match:
            if sync_success is None:
                status = "Awaiting Sync"
            elif sync_success:
                status = "Success"
            else:
                status = "Failed"
            destination_status_and_url[sd.name] = (status, remote_url)
        else:
            destination_status_and_url[sd.name] = (None, None)
    return destination_status_and_url
