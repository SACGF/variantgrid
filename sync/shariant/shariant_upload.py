from typing import Dict

from django.db.models.query_utils import Q
import json
import requests

from library.guardian_utils import admin_bot
from library.oauth import OAuthConnector
from library.utils import batch_iterator
from sync.models.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter
from classification.enums.classification_enums import ShareLevel
from classification.models import EvidenceKeyMap
from classification.models.classification import ClassificationModification
from classification.models.classification_utils import ClassificationJsonParams

SHARIANT_PRIVATE_FIELDS = ['patient_id', 'family_id', 'sample_id', 'patient_summary']

def insert_nones(data: Dict) -> Dict:
    for ekey in EvidenceKeyMap().all_keys:
        key = ekey.key
        if key not in data:
            data[key] = None
    return data


def sync_shariant_upload(sync_destination: SyncDestination, full_sync: bool = False) -> SyncRun:
    config = sync_destination.config
    shariant = OAuthConnector.shariant_oauth_connector()

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

        user = data.get('owner', {}).get('value')
        if not user:
            user = vcm.user.username

        user = user_mappings.get(user, user)
        data['owner'] = {'value': user}

        data = insert_nones(data)
        data = historical_converter.to_shariant(vcm, data)
        formatted_json['data'] = data
        formatted_json['publish'] = share_level

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
            q_not_somatic = ~Q(classification__evidence__allele_origin__value='Somatic')
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

        try:
            for batch in batch_iterator(qs, batch_size=10):

                json_to_send = {"records": [classification_to_json(vcm) for vcm in batch]}
                print(json.dumps(json_to_send))
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
