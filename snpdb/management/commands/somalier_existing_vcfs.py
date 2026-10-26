"""
Backfills somalier for every VCF without a successful extract: extracts fan out as celery tasks, then the new
extracts are grouped into ancestry batches of up to settings.SOMALIER["ancestry_batch_samples"] samples, then the
all-samples relate. Progress is read from the DB, so the command can be killed and re-run at any point and picks
up what is missing or failed (#1850). It fails if any extract or ancestry run errored, so the upgrader keeps the
step outstanding.
"""
import logging
import os
import time
from collections.abc import Callable, Iterable

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db.models import Count, Q, QuerySet
from django.utils import timezone

from library.vg.status import broker_queue_depths, task_queue_name
from snpdb.enums import ProcessingStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_somalier import (
    SomalierAllSamplesRelate,
    SomalierAncestryBatch,
    SomalierAncestryRun,
    SomalierConfig,
    SomalierRelatePairs,
    SomalierVCFExtract,
)
from snpdb.models.models_vcf import VCF
from snpdb.tasks.somalier_tasks import (
    STALE_PROCESSING_AGE,
    somalier_all_samples,
    somalier_ancestry_batch,
    somalier_vcf_extract,
)

POLL_SECONDS = 30
RUNNING_STATUSES = [ProcessingStatus.CREATED, ProcessingStatus.PROCESSING]


def batch_by_sample_count(extract_sample_counts: Iterable[tuple[int, int]], max_samples: int) -> list[list[int]]:
    """ (vcf_extract_id, num_samples) in order -> lists of extract ids of at most max_samples samples each,
        except that an extract bigger than that is a batch on its own """
    batches = []
    batch = []
    batch_samples = 0
    for vcf_extract_id, num_samples in extract_sample_counts:
        if batch and batch_samples + num_samples > max_samples:
            batches.append(batch)
            batch = []
            batch_samples = 0
        batch.append(vcf_extract_id)
        batch_samples += num_samples
    if batch:
        batches.append(batch)
    return batches


def _still_running_q(prefix: str) -> Q:
    """ A stage that's created or processing, and not left behind by a worker that died """
    return Q(**{f"{prefix}status__in": RUNNING_STATUSES,
                f"{prefix}modified__gt": timezone.now() - STALE_PROCESSING_AGE})


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--genome-build")
        parser.add_argument('--clear', action='store_true')
        parser.add_argument('--no-wait', action='store_true',
                            help="Only dispatch the extracts - re-run without it once the queue drains")

    def handle(self, *args, **options):
        if not settings.SOMALIER.get("enabled"):
            raise ValueError("settings.SOMALIER['enabled'] not enabled!")

        # Checked before --clear, as every extract would otherwise fail and the old ones are already gone
        somalier_bin = SomalierConfig().get_annotation("command")
        if not os.access(somalier_bin, os.X_OK):
            raise CommandError(f"somalier '{somalier_bin}' is missing or not executable - "
                               "run 'manage.py download_annotation_data'")

        if options.get("clear"):
            SomalierVCFExtract.objects.all().delete()
            SomalierAncestryBatch.objects.all().delete()  # These, and the relates, take their report dirs
            SomalierRelatePairs.objects.all().delete()
            SomalierAllSamplesRelate.objects.all().delete()

        vcf_qs = VCF.objects.all()
        if build_name := options.get("genome_build"):
            vcf_qs = vcf_qs.filter(genome_build=GenomeBuild.get_name_or_alias(build_name))

        needs_extract_q = Q(somaliervcfextract__isnull=True) | Q(somaliervcfextract__status=ProcessingStatus.ERROR)
        vcf_ids = list(vcf_qs.filter(needs_extract_q).order_by("pk").values_list("pk", flat=True))
        if not options.get("no_wait"):
            self._require_worker(somalier_vcf_extract)
            self._require_worker(somalier_ancestry_batch)
        logging.info("Dispatching somalier extract for %d VCFs", len(vcf_ids))
        for vcf_id in vcf_ids:
            somalier_vcf_extract.delay(vcf_id)

        if options.get("no_wait"):
            logging.info("Not waiting - run this command again without --no-wait once the queue has drained")
            return

        self._wait("extract", len(vcf_ids),
                   lambda: VCF.objects.filter(pk__in=vcf_ids).filter(
                       Q(somaliervcfextract__isnull=True) | _still_running_q("somaliervcfextract__")))

        SomalierAncestryRun.objects.filter(vcf_extract__vcf__in=vcf_qs, status=ProcessingStatus.ERROR).delete()
        extracts = SomalierVCFExtract.objects.filter(vcf__in=vcf_qs, status=ProcessingStatus.SUCCESS,
                                                     somalierancestryrun__isnull=True)
        extract_sample_counts = extracts.annotate(num_samples=Count("somaliersampleextract")) \
            .order_by("pk").values_list("pk", "num_samples")
        batches = batch_by_sample_count(extract_sample_counts, settings.SOMALIER["ancestry_batch_samples"])
        vcf_extract_ids = [pk for batch in batches for pk in batch]
        logging.info("Dispatching %d somalier ancestry batches for %d VCFs", len(batches), len(vcf_extract_ids))
        for batch in batches:
            somalier_ancestry_batch.delay(batch)

        self._wait("ancestry", len(vcf_extract_ids),
                   lambda: SomalierVCFExtract.objects.filter(pk__in=vcf_extract_ids).filter(
                       Q(somalierancestryrun__isnull=True) | _still_running_q("somalierancestryrun__")))

        somalier_all_samples.delay(force=True)
        logging.info("Dispatched the all samples relate")

        num_failed_extracts = SomalierVCFExtract.objects.filter(vcf__in=vcf_ids, status=ProcessingStatus.ERROR).count()
        num_failed_ancestry = SomalierAncestryRun.objects.filter(vcf_extract__in=vcf_extract_ids,
                                                                 status=ProcessingStatus.ERROR).count()
        if num_failed_extracts or num_failed_ancestry:
            raise CommandError(f"{num_failed_extracts} / {len(vcf_ids)} extracts and {num_failed_ancestry} / "
                               f"{len(vcf_extract_ids)} ancestry runs failed - see the Somalier stages on each "
                               "VCF's VCF Info tab, fix the cause and re-run (failed VCFs are retried)")

    @staticmethod
    def _require_worker(task):
        """ Waiting on a queue nobody consumes never ends, so stop before dispatching. --no-wait skips this
            as it exists to queue work ahead of a worker. """
        queue_name = task_queue_name(task.name)
        if broker_queue_depths([queue_name])[queue_name].get("consumers") == 0:
            raise CommandError(f"No celery worker is consuming queue '{queue_name}' ({task.name}) - start one "
                               f"(scripts/celery_start.sh on a dev box) and re-run, or use --no-wait")

    @staticmethod
    def _wait(stage: str, total: int, get_pending_qs: Callable[[], QuerySet]):
        """ get_pending_qs is called each poll, so a stage stuck since the last one becomes stale """
        while num_pending := get_pending_qs().count():
            logging.info("somalier %s: %d / %d done", stage, total - num_pending, total)
            time.sleep(POLL_SECONDS)
        logging.info("somalier %s: all %d done", stage, total)
