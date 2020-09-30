#!/usr/bin/env python3

import cyvcf2
from django.core.management.base import BaseCommand
import logging

from analysis.models.enums import NodeStatus
from analysis.vcf_files.bulk_vcf_count_inserter import BulkVCFCountInserter
from eventlog.models import create_event
from library.log_utils import get_traceback
from snpdb.models import VariantCollection


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('variant_collection_id', type=int)

    def handle(self, *args, **options):

        variant_collection_id = options['variant_collection_id']
        variant_collection = VariantCollection.objects.get(pk=variant_collection_id)

        logging.debug("Inserting variant_collection_id = %d", variant_collection_id)
        try:
            vcf_reader = cyvcf2.VCF("/dev/stdin")  # Must take a filename..
            bulk_inserter = BulkVCFCountInserter(variant_collection)

            for v in vcf_reader:
                bulk_inserter.process_entry(v)

            bulk_inserter.finish()  # Any leftovers
            variant_collection.count = bulk_inserter.rows_processed
            variant_collection.save()
        except Exception:
            details = get_traceback()
            logging.error(details)

            try:
                node = variant_collection.intersectioncache.node_version.node
                node.status = NodeStatus.ERROR
                errors = "Error inserting variants after bed intersection:\n"
                errors += details
                logging.error(errors)
                node.errors = errors
                node.save()
            except Exception as e:
                logging.error(e)
                create_event(None, name="stdin_to_variant_collection", details=details)
