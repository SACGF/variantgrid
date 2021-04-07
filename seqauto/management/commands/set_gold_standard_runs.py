#!/usr/bin/env python3
import re

from django.conf import settings
from django.core.management.base import BaseCommand
import logging
import os

from library.file_utils import file_to_array
from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_REGEX
from seqauto.models import EnrichmentKit, SequencingRun, GoldReference, ImportStatus
from seqauto.tasks.gold_summary_tasks import calculate_gold_summary


def get_enrichment_kit_gold_sequencing_runs(runs_for_enrichment_kit_qs):
    gold_qs = runs_for_enrichment_kit_qs.filter(gold_standard=True)
    return set(gold_qs.values_list("pk", flat=True))


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear', help="Clears all existing gold runs", action='store_true')
        parser.add_argument('--enrichment-kit', required=True)
        parser.add_argument('--enrichment-kit-version', type=int, default=1)
        parser.add_argument('gold_runs_filename')

    def handle(self, *args, **options):
        if options.get("clear"):
            logging.info("Clearing existing gold runs")
            SequencingRun.objects.all().update(gold_standard=False)

        gold_runs_filename = options["gold_runs_filename"]
        ek_name = options["enrichment_kit"]
        ek_version = options["enrichment_kit_version"]
        enrichment_kit = EnrichmentKit.objects.get(name=ek_name, version=ek_version)

        if not os.path.exists(gold_runs_filename):
            raise FileNotFoundError(gold_runs_filename)

        gold_runs = self._get_gold_runs(gold_runs_filename)
        runs_for_enrichment_kit_qs = SequencingRun.objects.filter(enrichment_kit=enrichment_kit)
        old_gold_for_enrichment_kit = get_enrichment_kit_gold_sequencing_runs(runs_for_enrichment_kit_qs)
        gold_runs_qs = SequencingRun.objects.all().filter(name__in=gold_runs)  # Get all in case was set to wrong kit
        logging.info("Setting enrichment_kit '%s' runs '%s' to this kit and gold", enrichment_kit, ','.join(gold_runs))
        # Also update enrichment_kit - as some runs weren't set properly in SampleSheet
        gold_runs_qs.update(gold_standard=True, enrichment_kit=enrichment_kit)
        gold_for_enrichment_kit = get_enrichment_kit_gold_sequencing_runs(runs_for_enrichment_kit_qs)

        missing = set(gold_runs) - gold_for_enrichment_kit
        if missing:
            logging.warning("*** Missing Gold runs ****")
            for sequencing_run_name in missing:
                logging.warning(sequencing_run_name)

        existing_gold_reference = GoldReference.objects.filter(enrichment_kit=enrichment_kit,
                                                               import_status=ImportStatus.SUCCESS)

        if old_gold_for_enrichment_kit != gold_for_enrichment_kit or not existing_gold_reference.exists():
            logging.info("Launched job 'calculate_gold_summary'")
            task = calculate_gold_summary.si(enrichment_kit.pk)  # @UndefinedVariable
            task.apply()

    @staticmethod
    def _get_gold_runs(gold_runs_filename):
        gold_runs = []
        for filename in file_to_array(gold_runs_filename):
            if filename.endswith("/"):
                filename = filename[:-1]
            gold_runs.append(os.path.basename(filename))
        return gold_runs
