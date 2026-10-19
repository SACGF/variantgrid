"""
    Inheritance mode stats for a Trio - the counts a TrioNode reads for denovo, autosomal recessive and
    dominant, X-linked recessive and compound het. They are written to the trio's cohort's CGC as
    aggregate rows under Trio.get_stats_filter_key, alongside the cohort's own stats, and computed from
    the same cohort rows (calculate_sample_stats), which know nothing about trios. A cohort recompute
    clears them, and the next TrioNode read misses and recomputes its trio (#2327).

    Entry points: calculate_trio_stats, enqueue_trio_stats_recompute (a TrioNode cache miss)
"""
import logging
from collections import defaultdict
from typing import Optional

import celery
from django.conf import settings
from django.db import connection, transaction

from annotation.models import AnnotationVersion
from annotation.tasks.calculate_sample_stats import (
    BOOLEAN_FLAGS_PARAMS,
    FLAG_BITS_SQL,
    SYMBOLIC_ALT_SQL,
    add_count,
    cohort_has_filters,
    cohort_stats_are_fresh,
    decode_flags,
    get_cohort_rows_queryset,
    get_sample_stats_code_version,
    new_counter,
    persist_cohort_stats,
)
from library.django_utils.database_utils import pg_settings
from snpdb.models import CohortGenotypeCollection, Zygosity
from snpdb.models.models_family import Trio


@celery.shared_task
def calculate_trio_stats_task(trio_id, annotation_version_id):
    trio = Trio.objects.get(pk=trio_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    calculate_trio_stats(trio, annotation_version)


def enqueue_trio_stats_recompute(trio: Trio, annotation_version: AnnotationVersion) -> None:
    """ Fire-and-forget recompute for a TrioNode that missed the cache. Idempotent -
        calculate_trio_stats early-exits if the rows are already fresh. """
    calculate_trio_stats_task.apply_async(
        args=[trio.pk, annotation_version.pk],
        queue="annotation_workers",
    )


def calculate_trio_stats(trio: Trio, annotation_version: AnnotationVersion):
    """ Count the trio's inheritance mode buckets and replace its rows on the cohort's CGC. Takes the same
        per-CGC lock as calculate_cohort_stats, and skips trios whose buckets are fresh - which includes a
        second trio with the same members and affected flags, as its keys are the same """
    code_version = get_sample_stats_code_version()
    filter_keys = trio.get_stats_filter_keys()
    cohort = trio.cohort
    with transaction.atomic():
        cgc = CohortGenotypeCollection.objects.select_for_update().get(pk=cohort.cohort_genotype_collection.pk)
        if cohort_stats_are_fresh(cgc, annotation_version, code_version, filter_keys, sample_ids=[]):
            logging.info("Trio stats already fresh for %s on %s; skipping recompute", trio, cgc)
            return
        has_filters = cohort_has_filters(cohort)
        aggregate_counters = _count_trio_stats(trio, cgc, annotation_version, has_filters)
        persist_cohort_stats(cgc, annotation_version, code_version, has_filters,
                             {}, aggregate_counters, filter_keys, samples=[], replace_all=False)


def _count_trio_stats(trio: Trio, cgc: CohortGenotypeCollection, annotation_version: AnnotationVersion,
                      has_filters: bool) -> dict[tuple[str, bool], list[int]]:
    """ (filter_key, passing_filter) → counter array for persist_cohort_stats. The cohort rows are grouped
        by the three members' zygosities in SQL, then Python applies the inheritance predicates and the
        compound het gene pairing (#2283) """
    aggregate_counters: dict[tuple[str, bool], list[int]] = defaultdict(new_counter)
    if (indexes := _get_packed_indexes(trio)) is None:
        return aggregate_counters

    def passing_filter_buckets(flags):
        return (False, True) if has_filters and flags.passing_filter else (False,)

    rows_sql, rows_params = get_cohort_rows_queryset(cgc, annotation_version, with_gene=True).query.sql_with_params()
    mum, dad, proband = (f"substr(v.samples_zygosity, {i + 1}, 1)" for i in indexes)  # 1-based
    with (pg_settings(local=True, jit="off", work_mem=settings.ANNOTATION_COHORT_STATS_WORK_MEM),
          connection.cursor() as cursor):
        # Only a proband HET can be compound het (_compound_het_source), so only those rows keep their gene
        cursor.execute(f"SELECT {mum}, {dad}, {proband}, CASE WHEN {proband} = %s THEN v.gene END, "
                       f"{SYMBOLIC_ALT_SQL}, {FLAG_BITS_SQL}, count(*) "
                       f"FROM ({rows_sql}) v GROUP BY 1, 2, 3, 4, 5, 6",
                       [Zygosity.HET, *BOOLEAN_FLAGS_PARAMS, *rows_params])
        _count_trio_modes(trio, cursor.fetchall(), aggregate_counters, passing_filter_buckets)
    return aggregate_counters


def _get_packed_indexes(trio: Trio) -> Optional[tuple[int, int, int]]:
    """ The trio's (mother, father, proband) packed indexes, or None when a member isn't in its cohort """
    packed_index_by_cohort_sample_id = {cs.pk: cs.cohort_genotype_packed_field_index
                                        for cs in trio.cohort.get_cohort_samples()}
    try:
        return (packed_index_by_cohort_sample_id[trio.mother_id],
                packed_index_by_cohort_sample_id[trio.father_id],
                packed_index_by_cohort_sample_id[trio.proband_id])
    except KeyError:
        return None


def _count_trio_modes(trio: Trio, trio_rows, aggregate_counters, passing_filter_buckets):
    """ trio_rows are (mum, dad, proband zygosity, gene_id, symbolic_alt, flag_bits, n). Each inheritance mode
        bucket counts like the cohort aggregate, classified by the proband's zygosity, without chrX fields.
        Compound het rows (proband HET from exactly one parent) only count for genes that have a
        mum-only and a dad-only source """
    key_by_mode = {mode: trio.get_stats_filter_key(mode) for mode in Trio.STATS_INHERITANCE_MODES}
    chet_key = key_by_mode["compound_het"]
    chet_gene_sources: dict[tuple[int, bool], set] = defaultdict(set)
    chet_pending = []
    for mum_zyg, dad_zyg, proband_zyg, gene_id, symbolic_alt, flag_bits, n in trio_rows:
        flags = decode_flags(symbolic_alt, flag_bits)
        modes = _trio_predicates(mum_zyg, dad_zyg, proband_zyg, flags.is_x,
                                 trio.mother_affected, trio.father_affected)
        chet_source = _compound_het_source(mum_zyg, dad_zyg, proband_zyg)
        for pf_pass in passing_filter_buckets(flags):
            for mode in modes:
                add_count(aggregate_counters[(key_by_mode[mode], pf_pass)], flags, proband_zyg, n, count_chrx=False)
            if chet_source is not None and gene_id is not None:
                chet_gene_sources[(gene_id, pf_pass)].add(chet_source)
                chet_pending.append((gene_id, pf_pass, flags, proband_zyg, n))

    for gene_id, pf_pass, flags, proband_zyg, n in chet_pending:
        if {"mum", "dad"} <= chet_gene_sources[(gene_id, pf_pass)]:
            add_count(aggregate_counters[(chet_key, pf_pass)], flags, proband_zyg, n, count_chrx=False)


def _trio_predicates(mum_zyg, dad_zyg, proband_zyg, is_x: bool, mother_affected, father_affected):
    """ Returns the set of inheritance modes the variant qualifies for under
        the trio's mom/dad/proband zygosities. Mirrors the predicates in
        analysis/models/nodes/sources/trio_node.py — keep in sync. """
    NO_VARIANT = (Zygosity.MISSING, Zygosity.HOM_REF)
    HAS_VARIANT = (Zygosity.HET, Zygosity.HOM_ALT)
    modes = set()
    # Denovo: mum & dad no variant; proband has variant
    if mum_zyg in NO_VARIANT and dad_zyg in NO_VARIANT and proband_zyg in HAS_VARIANT:
        modes.add("denovo")
    # Autosomal recessive: mum HET, dad HET, proband HOM_ALT
    if mum_zyg == Zygosity.HET and dad_zyg == Zygosity.HET and proband_zyg == Zygosity.HOM_ALT:
        modes.add("autosomal_recessive")
    # Autosomal dominant: each parent's expected zyg depends on affected status;
    # proband always has variant.
    mum_dominant = (mum_zyg in HAS_VARIANT) if mother_affected else (mum_zyg in NO_VARIANT)
    dad_dominant = (dad_zyg in HAS_VARIANT) if father_affected else (dad_zyg in NO_VARIANT)
    if mum_dominant and dad_dominant and proband_zyg in HAS_VARIANT:
        modes.add("autosomal_dominant")
    # X-linked recessive: mum HET, dad unconstrained, proband HOM_ALT, on chrX
    if is_x and mum_zyg == Zygosity.HET and proband_zyg == Zygosity.HOM_ALT:
        modes.add("x_linked_recessive")
    return modes


def _compound_het_source(mum_zyg, dad_zyg, proband_zyg) -> Optional[str]:
    """ The parent a proband HET is inherited from, when exactly one parent is HET """
    no_variant = (Zygosity.MISSING, Zygosity.HOM_REF)
    if proband_zyg == Zygosity.HET:
        if mum_zyg == Zygosity.HET and dad_zyg in no_variant:
            return "mum"
        if mum_zyg in no_variant and dad_zyg == Zygosity.HET:
            return "dad"
    return None
