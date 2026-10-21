"""
Backfills somalier for every VCF without an extract: extracts fan out as celery tasks, then the new extracts are
grouped into ancestry batches of up to settings.SOMALIER["ancestry_batch_samples"] samples, then the all-samples
relate. Progress is read from the DB, so the command can be killed and re-run at any point and picks up what
is missing (#1850).
"""
import logging
import time
from collections.abc import Callable, Iterable

from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import Count, Q, QuerySet
from django.utils import timezone

from snpdb.enums import ProcessingStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_somalier import (
    SomalierAllSamplesRelate,
    SomalierAncestryBatch,
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

        if options.get("clear"):
            SomalierVCFExtract.objects.all().delete()
            SomalierAncestryBatch.objects.all().delete()  # These, and the relates, take their report dirs
            SomalierRelatePairs.objects.all().delete()
            SomalierAllSamplesRelate.objects.all().delete()

        vcf_qs = VCF.objects.all()
        if build_name := options.get("genome_build"):
            vcf_qs = vcf_qs.filter(genome_build=GenomeBuild.get_name_or_alias(build_name))

        vcf_ids = list(vcf_qs.filter(somaliervcfextract__isnull=True).order_by("pk").values_list("pk", flat=True))
        logging.info("Dispatching somalier extract for %d VCFs", len(vcf_ids))
        for vcf_id in vcf_ids:
            somalier_vcf_extract.delay(vcf_id)

        if options.get("no_wait"):
            logging.info("Not waiting - run this command again without --no-wait once the queue has drained")
            return

        self._wait("extract", len(vcf_ids),
                   lambda: VCF.objects.filter(pk__in=vcf_ids).filter(
                       Q(somaliervcfextract__isnull=True) | _still_running_q("somaliervcfextract__")))

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

    @staticmethod
    def _wait(stage: str, total: int, get_pending_qs: Callable[[], QuerySet]):
        """ get_pending_qs is called each poll, so a stage stuck since the last one becomes stale """
        while num_pending := get_pending_qs().count():
            logging.info("somalier %s: %d / %d done", stage, total - num_pending, total)
            time.sleep(POLL_SECONDS)
        logging.info("somalier %s: all %d done", stage, total)
