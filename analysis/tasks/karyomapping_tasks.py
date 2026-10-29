from collections import Counter, defaultdict

import celery
from django.db import transaction
from django.utils import timezone

from analysis.models.models_karyomapping import (
    ContigKaryomappingCounts,
    GenomeKaryomappingCounts,
    KaryotypeBins,
)
from snpdb.enums import ImportStatus
from snpdb.models import Trio


def get_or_queue_genome_karyomapping_counts(trio: Trio) -> GenomeKaryomappingCounts:
    """ Counting reads the whole cohort, so a trio with none (eg deleted when binning changed) is counted when viewed,
        as is one whose count failed (eg the trio was made before its cohort genotypes) or never finished """
    genome_karyomapping_counts, queue = GenomeKaryomappingCounts.objects.get_or_create(trio=trio)
    if not queue and (retry_after := genome_karyomapping_counts.retry_after) and retry_after <= timezone.now():
        # Conditional on modified, so concurrent views requeue it once
        queue = GenomeKaryomappingCounts.objects.filter(
            pk=genome_karyomapping_counts.pk, modified=genome_karyomapping_counts.modified,
        ).update(import_status=ImportStatus.CREATED, modified=timezone.now())
        genome_karyomapping_counts.refresh_from_db()
    if queue:
        transaction.on_commit(create_genome_karyomapping_for_trio.si(trio.pk).apply_async)
    return genome_karyomapping_counts


@celery.shared_task
def create_genome_karyomapping_for_trio(trio_id):
    trio = Trio.objects.get(pk=trio_id)

    genome_karyomapping_counts, _ = GenomeKaryomappingCounts.objects.get_or_create(trio=trio)
    genome_karyomapping_counts.import_status = ImportStatus.IMPORTING
    genome_karyomapping_counts.save()
    try:
        _count_genome_karyomapping(genome_karyomapping_counts)
    except Exception:
        genome_karyomapping_counts.import_status = ImportStatus.ERROR
        genome_karyomapping_counts.save()
        raise


def _count_genome_karyomapping(genome_karyomapping_counts: GenomeKaryomappingCounts):
    trio = genome_karyomapping_counts.trio
    genome_karyomapping_counts.contigkaryomappingcounts_set.all().delete()

    contig_code_count = defaultdict(Counter)
    for variant_data, genotype_tuple in KaryotypeBins.iter_variant_and_genotypes(trio):
        if karyotype_bin := KaryotypeBins.get_bin(genotype_tuple):
            contig_id = variant_data[1]
            contig_code_count[contig_id][karyotype_bin] += 1

    genome_counts = Counter()
    contig_karyomapping_counts_list = []
    for contig_id, code_count in contig_code_count.items():
        genome_counts.update(code_count)

        kwargs = {"genome_karyomapping_counts": genome_karyomapping_counts,
                  "contig_id": contig_id}
        for code, count in code_count.items():
            kwargs[code.lower()] = count

        contig_karyomapping_counts_list.append(ContigKaryomappingCounts(**kwargs))

    ContigKaryomappingCounts.objects.bulk_create(contig_karyomapping_counts_list)

    for code in KaryotypeBins.COLLAPSED_BINS:
        setattr(genome_karyomapping_counts, code.lower(), genome_counts.get(code))
    genome_karyomapping_counts.import_status = ImportStatus.SUCCESS
    genome_karyomapping_counts.save()
