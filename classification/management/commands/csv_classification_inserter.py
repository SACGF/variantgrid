import logging
import time

import pandas as pd
from django.core.management import BaseCommand

from classification.enums import SubmissionSource
from classification.models import EvidenceKey
try:
    # Master code
    from classification.models.classification_inserter import BulkClassificationInserter
except ImportError:
    # VG3 SA Path Prod old code
    from classification.views.classification_view import BulkInserter as BulkClassificationInserter

from library.guardian_utils import admin_bot
from library.log_utils import log_traceback
from library.pandas_utils import df_nan_to_none
from library.utils import batch_iterator
from snpdb.models import Lab


class Command(BaseCommand):
    """
        Create classifications from ekey CSV files, @see https://github.com/SACGF/variantgrid/issues/1481

        TODO: (maybe optional and do from args)
            * set classification.sample (will be used in auto populate)
            * Auto populate - classification_auto_populate_fields
    """

    def add_arguments(self, parser):
        parser.add_argument('--max-records', type=int, default=1,  # TODO: Remove after testing ok
                            help="Limit number of records created")
        parser.add_argument('--lab', required=True, type=str,
                            help="Lab name")
        parser.add_argument('--variants-ekeys-csv', required=True,
                            help="Each row represents a new Classification, columns = Ekeys")
        parser.add_argument('--static-ekeys-csv', required=True,
                            help="These are applied to every classification")
        parser.add_argument('--sleep', type=int, default=0, required=False,
                            help="Sleep seconds between batches")



    def handle(self, *args, **options):
        max_records = options["max_records"]
        lab_name = options["lab"]
        variants_ekeys_csv = options["variants_ekeys_csv"]
        static_ekeys_csv = options["static_ekeys_csv"]
        sleep = options["sleep"]
        static_dict = self.get_static_keys(static_ekeys_csv)

        user = admin_bot()  # Or make it --user param?
        lab = Lab.objects.get(name=lab_name)

        records = self.iter_records(lab, static_dict, variants_ekeys_csv)

        count = 0
        first_batch = False
        for batch in batch_iterator(records, batch_size=50):
            inserter = BulkClassificationInserter(user=user, force_publish=True)
            try:
                # give the process 10 seconds to breath between batches of 50 classifications
                # in case we're downloading giant chunks of data
                if sleep and not first_batch:
                    logging.info("Sleeping %d secs" % sleep)
                    time.sleep(sleep)
                first_batch = False

                for record in batch:

                    inserter.insert(
                        record,
                        # record_id=record.pop("record_id"),
                        # submission_source=SubmissionSource.API,
                    )
                    count = count + 1
                    if count >= max_records:
                        break
            except Exception as e:
                log_traceback()
                raise
            finally:
                logging.info("Finish")
                inserter.finish()

            if count >= max_records:
                break

    @staticmethod
    def get_static_keys(static_ekeys_csv):
        df = pd.read_csv(static_ekeys_csv)
        row = df.loc[0]
        return row.to_dict()

    def get_internal_notes_ekey_name(self):
        # The ekey to store internal data differs between VG3 and master
        POTENTIAL_EKEYS = [
            "internal_use", # VG3 SA Path
            "review_comment",  # Master
        ]
        if ekey := EvidenceKey.objects.filter(pk__in=POTENTIAL_EKEYS).first():
            return ekey.key
        raise ValueError("No Ekey found")

    def iter_records(self, lab: Lab, static_dict, variants_ekeys_csv: str):
        # Are we going to do any validation? Code around the place uses:
        # known_keys = EvidenceKeyMap.instance()
        # known_keys.get(key).is_dummy
        # valid_evidence_keys = set(k.key for k in EvidenceKeyMap.cached().all_keys)

        internal_notes_ekey = self.get_internal_notes_ekey_name()
        df = pd.read_csv(variants_ekeys_csv)
        df = df_nan_to_none(df)

        for _, row in df.iterrows():
            data = static_dict.copy()
            data.update(row.to_dict())

            lab_record_id = data.pop("lab_record_id")
            id_str = f"{lab.group_name}/{lab_record_id}"  # Format is lab_id/]record_id[.version
            # Special case: column name differs from ekey name
            data[internal_notes_ekey] = {
                "value": data.pop("internal_use"),
                "note": data.pop("internal_use_note"),
            }
            for ekey in ["somatic:tumor_cellularity", "somatic:tmb_status", "somatic:msi_status"]:
                data[ekey] = {
                    "value": data.pop(ekey),
                    "note": data.pop(f"{ekey}_note"),
                }

            record = {
                "id": id_str,
                "upsert": data, # Whatever is left are popping
             }

            # logging.info(record)


            yield record
