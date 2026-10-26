"""
Recalculates ResolvedVariantInfo.p_hgvs / p_hgvs_source / p_hgvs_error (#1383), writing only those three columns.

Entry points:
- recalc_p_hgvs_for_pk_range (in process) and recalc_p_hgvs_for_pk_range_task, queued a batch at a time by the
  recalc_resolved_variant_info_p_hgvs command.
- recalc_p_hgvs_for_annotated_variants_task, queued when an annotation run completes (#2421): classifications are
  matched before VEP annotates a new variant, so their p.HGVS is refreshed once VEP's hgvs_p has landed.
"""
import itertools
import logging
from collections.abc import Iterator
from typing import Optional

import celery
from django.db.models import QuerySet

from classification.enums import PHGVSSource
from classification.models.classification_variant_info_models import (
    ResolvedVariantInfo,
    p_hgvs_converter,
)
from snpdb.models.models_genome import GenomeBuild

P_HGVS_FIELDS = ["p_hgvs", "p_hgvs_source", "p_hgvs_error"]
P_HGVS_BATCH_SIZE = 1000


def resolved_variant_infos_for_p_hgvs(genome_build_name: Optional[str] = None,
                                      only_missing: bool = False) -> QuerySet[ResolvedVariantInfo]:
    rvi_qs = ResolvedVariantInfo.objects.filter(transcript_version__isnull=False)
    if genome_build_name:
        rvi_qs = rvi_qs.filter(genome_build=GenomeBuild.get_name_or_alias(genome_build_name))
    if only_missing:
        rvi_qs = rvi_qs.filter(p_hgvs__isnull=True)
    return rvi_qs


def resolved_variant_infos_without_vep_p_hgvs(genome_build_name: str, min_variant_id: int,
                                              max_variant_id: int) -> QuerySet[ResolvedVariantInfo]:
    """ Rows on variants in the (inclusive) range whose p.HGVS didn't come from VEP """
    rvi_qs = resolved_variant_infos_for_p_hgvs(genome_build_name)
    rvi_qs = rvi_qs.filter(variant_id__gte=min_variant_id, variant_id__lte=max_variant_id)
    return rvi_qs.exclude(p_hgvs_source=PHGVSSource.VEP)


def pk_ranges(rvi_qs: QuerySet[ResolvedVariantInfo], batch_size: int) -> Iterator[tuple[int, int]]:
    """ Inclusive (min, max) pk ranges, each covering batch_size rows """
    pk_qs = rvi_qs.order_by("pk").values_list("pk", flat=True)
    for pks in itertools.batched(pk_qs.iterator(chunk_size=batch_size), batch_size):
        yield pks[0], pks[-1]


def recalc_p_hgvs_for_queryset(rvi_qs: QuerySet[ResolvedVariantInfo], retrieve_transcripts: bool = False) -> int:
    """ Returns the number of rows written """
    rvi_qs = rvi_qs.select_related("genome_build", "variant", "transcript_version",
                                   "allele_info__allele__clingen_allele")
    rvis = list(rvi_qs)
    for genome_build, build_rvis in itertools.groupby(sorted(rvis, key=lambda r: r.genome_build_id),
                                                      key=lambda r: r.genome_build):
        build_rvis = list(build_rvis)
        converter = p_hgvs_converter(genome_build, retrieve_transcripts=retrieve_transcripts)
        vep_hgvs_p_lookup = ResolvedVariantInfo.vep_hgvs_p_lookup_for(genome_build, build_rvis)
        for rvi in build_rvis:
            rvi.set_p_hgvs(converter=converter, vep_hgvs_p_lookup=vep_hgvs_p_lookup)
    ResolvedVariantInfo.objects.bulk_update(rvis, P_HGVS_FIELDS)
    return len(rvis)


def recalc_p_hgvs_for_pk_range(min_pk: int, max_pk: int, genome_build_name: Optional[str] = None,
                               only_missing: bool = False, retrieve_transcripts: bool = False) -> int:
    """ Returns the number of rows written """
    rvi_qs = resolved_variant_infos_for_p_hgvs(genome_build_name, only_missing).filter(pk__gte=min_pk, pk__lte=max_pk)
    return recalc_p_hgvs_for_queryset(rvi_qs, retrieve_transcripts=retrieve_transcripts)


@celery.shared_task(queue="db_workers")
def recalc_p_hgvs_for_pk_range_task(min_pk: int, max_pk: int, genome_build_name: Optional[str] = None,
                                    only_missing: bool = False, retrieve_transcripts: bool = False):
    num_processed = recalc_p_hgvs_for_pk_range(min_pk, max_pk, genome_build_name=genome_build_name,
                                               only_missing=only_missing, retrieve_transcripts=retrieve_transcripts)
    logging.info("ResolvedVariantInfo %d-%d: recalculated p.HGVS for %d rows", min_pk, max_pk, num_processed)


def recalc_p_hgvs_for_annotated_variants(genome_build_name: str, min_variant_id: int, max_variant_id: int) -> int:
    """ Returns the number of rows written """
    rvi_qs = resolved_variant_infos_without_vep_p_hgvs(genome_build_name, min_variant_id, max_variant_id)
    num_processed = 0
    for min_pk, max_pk in pk_ranges(rvi_qs, P_HGVS_BATCH_SIZE):
        num_processed += recalc_p_hgvs_for_queryset(rvi_qs.filter(pk__gte=min_pk, pk__lte=max_pk))
    return num_processed


@celery.shared_task(queue="db_workers")
def recalc_p_hgvs_for_annotated_variants_task(genome_build_name: str, min_variant_id: int, max_variant_id: int):
    num_processed = recalc_p_hgvs_for_annotated_variants(genome_build_name, min_variant_id, max_variant_id)
    logging.info("Variants %d-%d (%s) annotated: recalculated p.HGVS for %d ResolvedVariantInfo rows",
                 min_variant_id, max_variant_id, genome_build_name, num_processed)
