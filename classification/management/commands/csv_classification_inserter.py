import logging
import time

import pandas as pd
from django.core.management import BaseCommand

from classification.enums import SubmissionSource
from classification.models import EvidenceKeyMap, EvidenceKey
from classification.models.classification_inserter import BulkClassificationInserter
from library.guardian_utils import admin_bot
from library.pandas_utils import df_nan_to_none
from library.utils import batch_iterator


class Command(BaseCommand):
    """
        Create classifications from ekey CSV files, @see https://github.com/SACGF/variantgrid/issues/1481

        TODO:
            * "internal_use" Ekey to document metadata in VG3, but what to use in master? source_data

            * set classification.sample (will be used in auto populate)
            * Auto populate - classification_auto_populate_fields
    """

    def add_arguments(self, parser):
        parser.add_argument('--max-records', type=int, default=1,  # TODO: Remove after testing ok
                            help="Limit number of records created")
        parser.add_argument('--variants-ekeys-csv', required=True,
                            help="Each row represents a new Classification, columns = Ekeys")
        parser.add_argument('--static-ekeys-csv', required=True,
                            help="These are applied to every classification")

    def handle(self, *args, **options):
        max_records = options["max_records"]
        variants_ekeys_csv = options["variants_ekeys_csv"]
        static_ekeys_csv = options["static_ekeys_csv"]
        static_dict = self.get_static_keys(static_ekeys_csv)

        user = admin_bot()  # Or make it --user param?

        records = self.iter_records(static_dict, variants_ekeys_csv)

        count = 0
        first_batch = False
        for batch in batch_iterator(records, batch_size=50):
            inserter = BulkClassificationInserter(user=user, force_publish=True)
            try:
                # give the process 10 seconds to breath between batches of 50 classifications
                # in case we're downloading giant chunks of data
                if not first_batch:
                    time.sleep(10)
                first_batch = False

                for record in batch:
                    inserter.insert(
                        record,
                        submission_source=SubmissionSource.API,
                    )
                    count = count + 1
                    if count >= max_records:
                        break
            finally:
                logging.info("Finish")
                inserter.finish()

    @staticmethod
    def get_static_keys(static_ekeys_csv) -> dict[str, str]:
        df = pd.read_csv(static_ekeys_csv)
        row = df.loc[0]
        return row.to_dict()

    def get_internal_notes_ekey_name(self) -> str:
        # The ekey to store internal data differs between VG3 and master
        POTENTIAL_EKEYS = [
            "internal_use", # VG3 SA Path
            "review_comment",  # Master
        ]
        if ekey := EvidenceKey.objects.filter(pk__in=POTENTIAL_EKEYS).first():
            return ekey.key
        raise ValueError("No Ekey found")

    def iter_records(self, static_dict: dict[str,str], variants_ekeys_csv: str):
        internal_notes_ekey = self.get_internal_notes_ekey_name()
        known_keys = EvidenceKeyMap.instance()

        df = pd.read_csv(variants_ekeys_csv)
        df = df_nan_to_none(df)

        for _, row in df.iterrows():
            data = row.to_dict()
            lab_record_id = data.pop("lab_record_id")
            internal_use = data.pop("internal_use")
            internal_use_note = data.pop("internal_use_note")
            data[internal_notes_ekey] = {
                "value": internal_use,
                "note": internal_use_note,
            }

            record = {
                 "id": lab_record_id,
                 "upsert": data, # Whatever is left are popping
             }

            yield record


            # Are we going to do any validation? Code around the place uses:
            # known_keys.get(key).is_dummy
            # valid_evidence_keys = set(k.key for k in EvidenceKeyMap.cached().all_keys)
