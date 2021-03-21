#!/usr/bin/env python3
import re

from django.conf import settings
from django.core.management.base import BaseCommand
import logging
import os

from library.file_utils import file_to_array
from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_REGEX
from seqauto.models import EnrichmentKit, SequencingRun
from seqauto.tasks.gold_summary_tasks import calculate_gold_summary


def get_enrichment_kit_gold_sequencing_runs(runs_for_enrichment_kit_qs):
    gold_qs = runs_for_enrichment_kit_qs.filter(gold_standard=True)
    return set(gold_qs.values_list("pk", flat=True))


class Command(BaseCommand):
    TAU_NEW_HG38_KITS = ("idt_exome", "idt_gmp_focus")

    def handle(self, *args, **options):

        clear_existing_gold_standard_runs = True  # TODO: command line param?
        if clear_existing_gold_standard_runs:
            SequencingRun.objects.all().update(gold_standard=False)

        # For each enrichment_kit open QC file
        for enrichment_kit in EnrichmentKit.objects.all():
            if enrichment_kit.name == "roche_1k_disease":
                version_dir = f"version{enrichment_kit.version}"
                gold_runs_file_template = os.path.join(settings.SEQAUTO_GOLD_BASE_DIR, "%(enrichment_kit)s", "gold", version_dir, "gold_ref_lists_%(enrichment_kit)s.txt")
            elif enrichment_kit.name in self.TAU_NEW_HG38_KITS:
                version_dir = f"v2"
                gold_runs_file_template = os.path.join(settings.SEQAUTO_GOLD_BASE_DIR, "%(enrichment_kit)s", version_dir, "gold", "gold_run_lists.txt")
            else:
                gold_runs_file_template = os.path.join(settings.SEQAUTO_GOLD_BASE_DIR, "%(enrichment_kit)s", "gold", "gold_ref_lists_%(enrichment_kit)s.txt")

            gold_runs_filename = gold_runs_file_template % {"enrichment_kit": enrichment_kit.name}
            if os.path.exists(gold_runs_filename):
                gold_runs = self._get_gold_runs(enrichment_kit, gold_runs_filename)
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


    @staticmethod
    def _get_gold_runs(enrichment_kit, gold_runs_filename):
        """ Need to deal with how this was stored on /tau over time """
        gold_runs = file_to_array(gold_runs_filename)
        if enrichment_kit.name in Command.TAU_NEW_HG38_KITS:
            fixed_runs = []
            for filename in gold_runs:
                if filename.endswith("/"):
                    filename = filename[:-1]
                basename = os.path.basename(filename)
                if m := re.match(f".*({SEQUENCING_RUN_REGEX})", basename):
                    sequencing_run_name = m.group(1)
                    sequencing_run_name += f"_{enrichment_kit.name}
                    fixed_runs.append(sequencing_run_name)
                else:
                    raise ValueError("Could not extract sequencing run out of {basename}")

        return gold_runs