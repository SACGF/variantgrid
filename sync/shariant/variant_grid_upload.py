import socket
from typing import Iterable, TypeVar, Union

from django.db.models import QuerySet

from classification.enums.classification_enums import ShareLevel
from classification.models import EvidenceKeyMap
from classification.models.classification import ClassificationModification
from classification.models.classification_utils import ClassificationJsonParams
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from library.oauth import ServerAuth
from sync.models.models import SyncDestination
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter
from sync.shariant.query_json_filter import QueryJsonFilter
from sync.sync_runner import register_sync_runner, SyncRunner, SyncRunInstance

# add variant_type to private fields as the key has been deprecated
SHARIANT_PRIVATE_FIELDS = [
    'age_units', 'dob', 'family_id', 'internal_use', 'patient_id', 'patient_summary', 'sample_id', 'variant_type'
]


def insert_nones(data: dict) -> dict:
    for e_key in EvidenceKeyMap.instance().all_keys:
        key = e_key.key
        if key not in data:
            data[key] = None
    return data


T = TypeVar("T")


def batch_iterator_end(iterable: Iterable[T], batch_size: int = 10) -> Iterable[Union[list[T], bool]]:
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


@register_sync_runner(config={"type": {"shariant", "variantgrid"}, "direction": "upload"})
class VariantGridUploadSyncer(SyncRunner):

    def __init__(self, sync_destination: SyncDestination):
        self.sync_destination = sync_destination

        config = self.sync_destination.config
        self.filters = config.get('filters', {})
        mapping = config.get('mapping', {})

        self.shariant = ServerAuth.for_sync_details(sync_destination.sync_details)
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
            qs = qs.filter(q)

        if not full_sync:
            qs = ClassificationModificationSyncRecord.filter_out_synced(
                qs=qs,
                destination=self.sync_destination
            )

        return qs

    def classification_to_json(self, vcm: ClassificationModification) -> dict:
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

    def sync(self, sync_run_instance: SyncRunInstance):
        qs = self.records_to_sync(full_sync=sync_run_instance.full_sync)

        rows_uploaded = 0

        if not qs.exists():
            sync_run_instance.run_completed(had_records=False)
        else:
            if max_rows := sync_run_instance.max_rows:
                qs = qs[:max_rows]

            site_name = socket.gethostname().lower().split('.')[0].replace('-', '')
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
                other_variant_grid = sync_run_instance.server_auth()

                response = other_variant_grid.post(
                    url_suffix='classification/api/classifications/v2/record/',
                    json=json_to_send,
                    timeout=MINUTE_SECS,
                )
                response.raise_for_status()
                # results are sent back in an array in the same order they were sent up
                results = response.json().get('results')

                for record, result in zip(batch, results):
                    rows_uploaded += 1
                    ClassificationModificationSyncRecord.objects.create(
                        run=sync_run_instance.sync_run,
                        classification_modification=record,
                        meta=result
                    )
            sync_run_instance.run_completed(had_records=True)
