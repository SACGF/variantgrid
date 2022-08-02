import socket
from typing import Dict, Optional, Iterable, List, TypeVar, Union

import requests
from django.db.models.query_utils import Q

from classification.enums.classification_enums import ShareLevel
from classification.models import EvidenceKeyMap
from classification.models.classification import ClassificationModification
from classification.models.classification_utils import ClassificationJsonParams
from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector
from sync.models.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter

# add variant_type to private fields as the key has been deprecated
SHARIANT_PRIVATE_FIELDS = ['patient_id', 'family_id', 'sample_id', 'patient_summary', 'internal_use', 'variant_type', 'age_units']


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
    :param iterable: Iteratble of T
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


def sync_shariant_upload(sync_destination: SyncDestination, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
    config = sync_destination.config
    shariant = OAuthConnector.shariant_oauth_connector(config["sync_details"])
    filters = config.get('filters', {})
    mapping = config.get('mapping', {})
    lab_mappings = mapping.get('labs', {})
    share_level_mappings = mapping.get('share_levels', {})
    user_mappings = mapping.get('users', {})
    historical_converter = HistoricalEKeyConverter()

    def classification_to_json(vcm: ClassificationModification) -> Dict:
        raw_json = vcm.as_json(ClassificationJsonParams(current_user=admin_bot(), include_data=True, include_messages=False, strip_complicated=True, api_version=2))
        formatted_json = {}

        lab = vcm.classification.lab
        mapped_lab_name = lab_mappings.get(lab.group_name)
        if mapped_lab_name is True:
            mapped_lab_name = lab.group_name

        share_level = vcm.share_level
        share_level = share_level_mappings.get(share_level, share_level)

        # might need to map the the lab group names if we don't map them all into shariant
        formatted_json['id'] = mapped_lab_name + '/' + vcm.classification.lab_record_id
        data = raw_json.get('data')
        for dont_share in SHARIANT_PRIVATE_FIELDS:
            data.pop(dont_share, None)

        # no need to screw around with owner
        user = data.get('owner', {}).get('value')
        if not user:
            user = vcm.user.username

        user = user_mappings.get(user, user)
        data['owner'] = {'value': user}

        data = insert_nones(data)
        data = historical_converter.to_shariant(vcm, data)
        formatted_json['overwrite'] = data
        formatted_json['publish'] = share_level
        # only allow for withdrawing, doesn't allow for un-withdrawing
        # but Shariant should get a notification that a withdrawn classification is being re-sent
        if raw_json.get('withdrawn', False):
            formatted_json['delete'] = True

        return formatted_json

    # note that more mapping is done during the upload call

    qs = ClassificationModification.objects.filter(
        is_last_published=True,
        share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
        classification__lab__group_name__in=lab_mappings.keys()
    )
    # Originally I used filter/exclude dicts which were turned into queryset kwargs allowing
    # Arbitrary filtering, but this doesn't work on JSONB fields where
    # qs.exclude(classification__evidence__allele_origin__value='Somatic') returns 0 results
    # Use 1 off special case code to handle this, with config like 'filters' : {"somatic" : false}
    # If you end up making more filters, perhaps consider adopting arbitrary filter construction
    # possibly by pulling code out of JQGrid
    if filters:
        somatic = filters.get("somatic", True)
        if not somatic:
            q_allele_origin_null = Q(classification__evidence__allele_origin__isnull=True)
            q_not_somatic = ~Q(classification__evidence__allele_origin__value='somatic')
            qs = qs.filter(q_allele_origin_null | q_not_somatic)

    if not full_sync:
        qs = ClassificationModificationSyncRecord.filter_out_synced(
            qs=qs,
            destination=sync_destination
        )

    """
    Provide a QuerySet of ClassificationModifications
    """
    rows_uploaded = 0
    run = SyncRun(destination=sync_destination, status=SyncStatus.IN_PROGRESS)
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
                    "records": [classification_to_json(vcm) for vcm in batch],
                    "import_id": site_name
                }
                if finished:
                    json_to_send["status"] = "complete"
                # print(json.dumps(json_to_send))
                auth = shariant.auth()

                response = requests.post(
                    shariant.url('classification/api/classifications/v2/record/'),
                    auth=auth,
                    json=json_to_send
                )
                response.raise_for_status()
                # results are sent back in an array in the same order they were sent up
                results = response.json().get('results')

                for record, result in zip(batch, results):
                    rows_uploaded += 1
                    ClassificationModificationSyncRecord.objects.create(run=run, classification_modification=record, meta=result)

            run.status = SyncStatus.SUCCESS
            run.save()
        finally:

            # bit of a hack so I can let the exception just fall through to rollbar
            # while still updating the status of the run
            if run.status == SyncStatus.IN_PROGRESS:
                run.status = SyncStatus.FAILED
                run.save()

    return run
