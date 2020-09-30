from django.conf import settings
import logging

from library.file_utils import IteratorFile
from upload.vcf.sql_copy_files import sql_copy_csv_file


class BulkVCFCountInserter:
    COLUMNS = ['variant_id', 'variant_collection_id']

    def __init__(self, variant_collection):
        self.variant_collection = variant_collection
        self.variant_collection_id_str = str(variant_collection.pk)
        self.table_name = variant_collection.get_partition_table()
        self.rows_processed = 0
        self.variant_collection_record_list = []

    def process_entry(self, v):
        if len(self.variant_collection_record_list) >= settings.SQL_BATCH_INSERT_SIZE:
            self.bulk_insert()

        variant_id = str(v.INFO["variant_id"])
        line = ','.join((variant_id, self.variant_collection_id_str)) + '\n'
        self.variant_collection_record_list.append(line)
        self.rows_processed += 1

    def finish(self):
        self.bulk_insert()

    def bulk_insert(self):
        num_records = len(self.variant_collection_record_list)
        logging.info("bulk_insert %d records into %s", num_records, self.table_name)

        if num_records:

            def csv_iterator():
                for line in self.variant_collection_record_list:
                    yield line

            f = IteratorFile(csv_iterator())

            sql_copy_csv_file(f, self.table_name, BulkVCFCountInserter.COLUMNS)
            self.variant_collection_record_list = []
