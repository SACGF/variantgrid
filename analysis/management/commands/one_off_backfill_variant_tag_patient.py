"""
Fill in VariantTag.patient for taggings made before the field existed.

The patient is who a tagging is about where the node knows the person but not which of their samples - a
patient / specimen / extraction level source node whose arms come from several callers. It is what puts such a
tagging into the patient's classify queue and narrows the classify form's sample dropdown, so older taggings
need the same answer the tag-time path now records.

A tagging that already has a sample takes that sample's patient, since the sample is the more specific answer
and the node walk can be ambiguous where the sample was set from a wizard node. Anything still unknown is
left null and falls back to the analysis' single patient at display time.

Taggings are done an analysis at a time so the node graph - which is what the answer actually comes from - is
loaded once for the thousands of taggings that share it. Only taggings without a patient are selected, so an
interrupted run picks up where it stopped.

@see https://github.com/SACGF/variantgrid/issues/1854
"""
import logging
import time
from datetime import timedelta

from django.core.management import BaseCommand

from analysis.models import Analysis, VariantTag
from analysis.variant_tag_operations import get_proband_by_node_id
from patients.models_enums import SampleSourceLevel
from patients.sample_grouping import get_patient_for_source

BATCH_SIZE = 1000
PROGRESS_EVERY = 10_000  # Taggings between progress lines - analyses vary from a handful to thousands


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would be set without changing anything")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        base_qs = VariantTag.objects.filter(patient__isnull=True, analysis__isnull=False, node__isnull=False)
        analysis_ids = list(base_qs.order_by("analysis_id").values_list("analysis_id", flat=True).distinct())
        total = base_qs.count()  # Taggings still without a patient - a resumed run counts what is left to do
        self.stdout.write(f"{total} taggings without a patient, across {len(analysis_ids)} analyses")
        if not total:
            return

        start = time.time()
        filled = 0
        ambiguous = 0
        last_logged = 0
        for i, analysis_id in enumerate(analysis_ids, start=1):
            analysis = Analysis.objects.get(pk=analysis_id)
            proband_by_node_id = get_proband_by_node_id(analysis)

            to_update = []
            qs = base_qs.filter(analysis=analysis).select_related("sample__extraction__specimen")
            for variant_tag in qs.order_by("pk").iterator():
                patient = self._patient_for_tag(variant_tag, proband_by_node_id)
                if patient is None:
                    ambiguous += 1
                    continue
                variant_tag.patient = patient
                to_update.append(variant_tag)
                if len(to_update) >= BATCH_SIZE:
                    filled += self._write(to_update, dry_run)
                    to_update = []
            filled += self._write(to_update, dry_run)

            if filled + ambiguous - last_logged >= PROGRESS_EVERY or i == len(analysis_ids):
                last_logged = filled + ambiguous
                self._log_progress(analysis_id, i, len(analysis_ids), filled, ambiguous, total, start)

        prefix = "Would set" if dry_run else "Set"
        self.stdout.write(f"{prefix} patient on {filled} taggings, left {ambiguous} ambiguous")

    @staticmethod
    def _patient_for_tag(variant_tag: VariantTag, proband_by_node_id: dict):
        """ The tagging's own sample first - it says which person more precisely than a re-walk of the graph,
            which has moved on since the tag was made """
        if variant_tag.sample_id:
            if patient := get_patient_for_source(SampleSourceLevel.SAMPLE, variant_tag.sample):
                return patient
        if proband := proband_by_node_id.get(variant_tag.node_id):
            return proband.patient
        return None

    @staticmethod
    def _log_progress(analysis_id: int, analyses_done: int, analyses_total: int,
                      filled: int, ambiguous: int, total: int, start: float):
        done = filled + ambiguous
        elapsed = time.time() - start
        rate = done / elapsed if elapsed else 0
        remaining = timedelta(seconds=int((total - done) / rate)) if rate else "?"
        logging.info("analysis %s (%d/%d) - %d/%d taggings (%.1f%%), %d filled / %d ambiguous, "
                     "%.0f/sec, ~%s left",
                     analysis_id, analyses_done, analyses_total, done, total, 100 * done / total,
                     filled, ambiguous, rate, remaining)

    @staticmethod
    def _write(variant_tags: list[VariantTag], dry_run: bool) -> int:
        if variant_tags and not dry_run:
            VariantTag.objects.bulk_update(variant_tags, ["patient"])
        return len(variant_tags)
