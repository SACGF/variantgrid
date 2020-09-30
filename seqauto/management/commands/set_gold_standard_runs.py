#!/usr/bin/env python3

from django.conf import settings
from django.core.management.base import BaseCommand
import logging
import os

from library.file_utils import file_to_array
from seqauto.models import EnrichmentKit, SequencingRun
from seqauto.tasks.gold_summary_tasks import calculate_gold_summary


def get_enrichment_kit_gold_sequencing_runs(runs_for_enrichment_kit_qs):
    gold_qs = runs_for_enrichment_kit_qs.filter(gold_standard=True)
    return set(gold_qs.values_list("pk", flat=True))


class Command(BaseCommand):

    def handle(self, *args, **options):

        clear_existing_gold_standard_runs = True  # TODO: command line param?
        if clear_existing_gold_standard_runs:
            SequencingRun.objects.all().update(gold_standard=False)

        # For each enrichment_kit open QC file
        for enrichment_kit in EnrichmentKit.objects.all():
            if enrichment_kit.name == "roche_1k_disease":
                version_dir = f"version{enrichment_kit.version}"
                gold_runs_file_template = os.path.join(settings.SEQAUTO_QC_BASE_DIR, "%(enrichment_kit)s", "gold", version_dir, "gold_ref_lists_%(enrichment_kit)s.txt")
            else:
                gold_runs_file_template = os.path.join(settings.SEQAUTO_QC_BASE_DIR, "%(enrichment_kit)s", "gold", "gold_ref_lists_%(enrichment_kit)s.txt")
            gold_runs_filename = gold_runs_file_template % {"enrichment_kit": enrichment_kit.name}
            if os.path.exists(gold_runs_filename):
                gold_runs = file_to_array(gold_runs_filename)
                runs_for_enrichment_kit_qs = SequencingRun.objects.filter(samplesheet__sequencingsample__enrichment_kit=enrichment_kit).distinct()
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

                if old_gold_for_enrichment_kit != gold_for_enrichment_kit:
                    logging.info("Launched job 'calculate_gold_summary'")
                    task = calculate_gold_summary.si(enrichment_kit.pk)  # @UndefinedVariable
                    task.apply_async()
            else:
                logging.warning(f"Couldn't open gold runs filename: {gold_runs_filename}")
