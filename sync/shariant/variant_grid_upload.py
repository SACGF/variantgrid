import copy
import logging
import socket
import time

import requests
from collections.abc import Iterable
from typing import Any, Optional, TypeVar, Union
from urllib.parse import quote, urljoin

from django.db.models import QuerySet
from django.urls import reverse

from classification.enums.classification_enums import ShareLevel, SpecialEKeys
from classification.models.classification import ClassificationModification
from classification.models.classification_utils import ClassificationJsonParams
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from library.constants import MINUTE_SECS
from library.guardian_utils import admin_bot
from snpdb.models import Lab
from sync.models.models import SyncDestination
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter
from sync.shariant.query_json_filter import SOMATIC_FILTER_KEY, QueryJsonFilter
from sync.sync_runner import ClassificationUploadSyncRunner, SyncRunInstance, register_sync_runner

# add variant_type to private fields as the key has been deprecated
SHARIANT_PRIVATE_FIELDS = [
    'age_units', 'dob', 'family_id', 'internal_use', 'patient_id', 'patient_summary', 'sample_id', 'variant_type'
]

# server-side processing of a 50 record batch can exceed the default minute
UPLOAD_TIMEOUT_SECS = 5 * MINUTE_SECS
UPLOAD_ATTEMPTS = 3
UPLOAD_RETRY_DELAY_SECS = 30


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
class VariantGridUploadSyncer(ClassificationUploadSyncRunner):

    def __init__(self):
        self.sync_destination: Optional[SyncDestination] = None
        self.filters = {}
        self.filter_labels = {}
        self.remote_lab_record_url = False
        self.lab_mappings = {}
        self.share_level_mappings = {}
        self.user_mappings = {}
        self.historical_converter = HistoricalEKeyConverter()

    def configure(self, sync_destination: SyncDestination):
        self.sync_destination = sync_destination

        config = sync_destination.config
        self.filters = config.get('filters', {})
        self.filter_labels = config.get('filter_labels', {})
        # only true once the remote has been upgraded to a version serving view_classification_lab_record
        self.remote_lab_record_url = config.get('remote_lab_record_url', False)
        mapping = config.get('mapping', {})

        self.lab_mappings = mapping.get('labs', {})
        self.share_level_mappings = mapping.get('share_levels', {})
        self.user_mappings = mapping.get('users', {})

    def mapped_lab_name(self, lab: Lab) -> Optional[str]:
        mapped_lab_name = self.lab_mappings.get(lab.group_name)
        if mapped_lab_name is True:
            mapped_lab_name = lab.group_name
        return mapped_lab_name or None

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

    def exclusion_reasons(self, cm: ClassificationModification) -> list[str]:
        reasons = []
        if cm.share_level not in ShareLevel.DISCORDANT_LEVEL_KEYS:
            share_level_label = ShareLevel(cm.share_level).label
            reasons.append(f"Shared with {share_level_label} only, uploads once shared more widely")

        lab = cm.classification.lab
        if not self.mapped_lab_name(lab):
            reasons.append(f"{lab.name} is not configured to share with {self.sync_destination}")

        if self.filters:
            single_record_qs = ClassificationModification.objects.filter(pk=cm.pk)
            json_filter = QueryJsonFilter.classification_value_filter()
            for key, value in self.filters.items():
                # convert_to_q_w_key strips nulls out of list values, so give it a copy to chew on
                q = json_filter.convert_to_q_w_key(key, copy.deepcopy(value))
                if not single_record_qs.filter(q).exists():
                    reasons.append(self.describe_filter_clause(key, value, cm))

        return reasons

    def describe_filter_clause(self, key: str, value: Any, cm: ClassificationModification) -> str:
        if label := self.filter_labels.get(key):
            return label

        e_key_map = EvidenceKeyMap.cached()
        if key == SOMATIC_FILTER_KEY:
            # the hardcoded "somatic" clause is an allele origin test in disguise
            e_key = e_key_map.get(SpecialEKeys.ALLELE_ORIGIN)
            return f"{self._describe_current_value(e_key, cm)}, shares unless set to {e_key.pretty_value('somatic')}"

        if key in e_key_map:
            e_key = e_key_map.get(key)
            description = self._describe_current_value(e_key, cm)
            if permitted := self._describe_permitted_values(e_key, value):
                description += f", {permitted}"
            return description

        return f"{EvidenceKey.pretty_label_from_string(key)} excludes this record"

    @staticmethod
    def _describe_current_value(e_key: EvidenceKey, cm: ClassificationModification) -> str:
        return f"{e_key.pretty_label} is {e_key.pretty_value(cm.get(e_key.key)) or 'blank'}"

    @staticmethod
    def _describe_permitted_values(e_key: EvidenceKey, value: Any) -> Optional[str]:
        if isinstance(value, list):
            values = list(value)
            excludes = bool(values) and values[0] == '!'
            if excludes:
                values = values[1:]
            pretty = ", ".join(e_key.pretty_value(v) or 'blank' for v in values)
            return f"shares unless set to {pretty}" if excludes else f"shares when set to {pretty}"
        if isinstance(value, str):
            return f"shares when set to {e_key.pretty_value(value) or 'blank'}"
        return None

    def remote_url_for(self, cm: ClassificationModification) -> Optional[str]:
        mapped_lab_name = self.mapped_lab_name(cm.classification.lab)
        if not mapped_lab_name:
            return None

        lab_record_id = cm.classification.lab_record_id
        if self.remote_lab_record_url:
            path = reverse('view_classification_lab_record',
                           kwargs={'classification_ref': f"{mapped_lab_name}/{lab_record_id}"})
        else:
            # remote is on a version without view_classification_lab_record, so search for it instead
            path = f"variantopedia/search?search={quote(lab_record_id)}"
        return urljoin(self.sync_destination.sync_details["host"], path)

    def classification_to_json(self, vcm: ClassificationModification) -> dict:
        raw_json = vcm.as_json(ClassificationJsonParams(current_user=admin_bot(), include_data=True, include_messages=False, strip_complicated=True, api_version=2))
        formatted_json = {}

        lab = vcm.classification.lab
        mapped_lab_name = self.mapped_lab_name(lab)

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
        self.configure(sync_run_instance.sync_destination)
        qs = self.records_to_sync(full_sync=sync_run_instance.full_sync)

        rows_uploaded = 0
        record_count = qs.count()
        logging.info("%s: %d record(s) to upload", sync_run_instance.sync_destination, record_count)

        if not record_count:
            sync_run_instance.run_completed(had_records=False)
        else:
            if max_rows := sync_run_instance.max_rows:
                qs = qs[:max_rows]
                record_count = min(record_count, max_rows)

            site_name = socket.gethostname().lower().split('.')[0].replace('-', '')
            for batch, finished in batch_iterator_end(qs, batch_size=50):
                # providing import_id and status:complete when we're done lets Shariant know
                # when the upload has been completed
                json_start = time.time()
                json_to_send = {
                    "records": [self.classification_to_json(vcm) for vcm in batch],
                    "import_id": site_name
                }
                if finished:
                    json_to_send["status"] = "complete"
                json_duration = time.time() - json_start
                other_variant_grid = sync_run_instance.server_auth()

                post_start = time.time()
                # records are upserts keyed by lab/lab_record_id, so re-sending a batch after a
                # timeout/connection drop is safe even if the server processed the first attempt
                for attempt in range(UPLOAD_ATTEMPTS):
                    try:
                        response = other_variant_grid.post(
                            url_suffix='classification/api/classifications/v2/record/',
                            json=json_to_send,
                            timeout=UPLOAD_TIMEOUT_SECS,
                        )
                        response.raise_for_status()
                        break
                    except (requests.ConnectionError, requests.Timeout) as e:
                        retries_left = UPLOAD_ATTEMPTS - attempt - 1
                        if not retries_left:
                            raise
                        logging.warning("%s: %r - %d retries left", sync_run_instance.sync_destination, e, retries_left)
                        time.sleep(UPLOAD_RETRY_DELAY_SECS)
                post_duration = time.time() - post_start
                # results are sent back in an array in the same order they were sent up
                results = response.json().get('results')

                for record, result in zip(batch, results):
                    rows_uploaded += 1
                    ClassificationModificationSyncRecord.objects.create(
                        run=sync_run_instance.sync_run,
                        classification_modification=record,
                        meta=result
                    )
                logging.info("%s: uploaded %d/%d records (batch of %d: JSON build %.1fs, POST %.1fs)",
                             sync_run_instance.sync_destination, rows_uploaded, record_count,
                             len(batch), json_duration, post_duration)
            sync_run_instance.run_completed(had_records=True)
