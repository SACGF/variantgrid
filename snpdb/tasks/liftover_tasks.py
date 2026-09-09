import itertools
import logging
from collections.abc import Iterable

import celery
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet

from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Allele, GenomeBuild, ImportSource


@celery.shared_task
def liftover_alleles(username, genome_build_name, retry_conversion_tool: str = None):
    """ Queues a task per batch of alleles - the db_workers pool caps how many run at once

        retry_conversion_tool (an AlleleConversionTool value) restricts this to the alleles that tool has
        failed on, and re-attempts it for them - @see create_liftover_pipelines """
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    allele_qs = _alleles_to_liftover(genome_build, retry_conversion_tool)

    for other_build in GenomeBuild.builds_with_annotation():
        if not GenomeBuild.is_equivalent(genome_build, other_build):
            num_batches = 0
            for min_allele_id, max_allele_id in _allele_id_batches(allele_qs):
                liftover_allele_batch.si(username, genome_build.name, other_build.name,
                                         min_allele_id, max_allele_id, retry_conversion_tool).apply_async()
                num_batches += 1
            logging.info("Queued %d liftover batches from %s to %s", num_batches, other_build, genome_build)


def _alleles_to_liftover(genome_build: GenomeBuild, retry_conversion_tool: str = None) -> QuerySet[Allele]:
    if retry_conversion_tool:
        return Allele.failed_liftover_for_build(genome_build, retry_conversion_tool)
    return Allele.missing_variants_for_build(genome_build)


def _allele_id_batches(allele_qs: QuerySet) -> Iterable[tuple[int, int]]:
    """ Inclusive (min, max) allele id ranges, each covering LIFTOVER_BATCH_SIZE alleles """
    batch_size = settings.LIFTOVER_BATCH_SIZE
    allele_id_qs = allele_qs.order_by("pk").values_list("pk", flat=True)
    for allele_ids in itertools.batched(allele_id_qs.iterator(chunk_size=batch_size), batch_size):
        yield allele_ids[0], allele_ids[-1]


@celery.shared_task
def liftover_allele_batch(username, genome_build_name, other_build_name, min_allele_id, max_allele_id,
                          retry_conversion_tool: str = None):
    user = User.objects.get(username=username)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    other_build = GenomeBuild.get_name_or_alias(other_build_name)
    # Re-query, so any alleles lifted over since the batches were worked out drop out
    alleles = _alleles_to_liftover(genome_build, retry_conversion_tool).filter(pk__gte=min_allele_id,
                                                                              pk__lte=max_allele_id)
    retry_conversion_tools = [retry_conversion_tool] if retry_conversion_tool else []
    logging.info("creating liftover pipelines from %s to %s for alleles %d-%d",
                 other_build, genome_build, min_allele_id, max_allele_id)
    create_liftover_pipelines(user, alleles, ImportSource.WEB, other_build, [genome_build],
                              retry_conversion_tools=retry_conversion_tools)
