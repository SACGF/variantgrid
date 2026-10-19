import logging
import time
from collections import defaultdict
from functools import cache
from typing import NamedTuple, Optional

import celery
import numpy as np
from django.conf import settings
from django.db import connection, transaction
from django.db.models import F
from django.db.models.functions import Length

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import (
    AnnotationVersion,
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models_annotation_run import VCFAnnotationStats
from eventlog.models import create_event
from library.django_utils import thread_safe_unique_together_get_or_create
from library.django_utils.database_utils import pg_settings
from library.enums.log_level import LogLevel
from library.genomics.vcf_enums import GeneLevelSymbolicAlt, VariantClass, VCFSymbolicAllele
from library.log_utils import get_traceback
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    ImportStatus,
    SampleStatsCodeVersion,
    Sequence,
    Variant,
    VCFLengthStats,
    VCFLengthStatsCollection,
    Zygosity,
)
from snpdb.models.models import CodeVersion
from snpdb.models.models_family import Trio

# Bump when a change alters what the stats count - readers treat rows from older versions as stale
# and enqueue a recompute (@see snpdb/views/views_data.py:_sample_stats)
SAMPLE_STATS_CODE_VERSION = 5


@celery.shared_task
def calculate_vcf_stats(vcf_id, annotation_version_id):
    vcf = VCF.objects.get(pk=vcf_id)
    try:
        annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
        msg = f"Build for annotation: {annotation_version.genome_build} and VCF: {vcf.genome_build} match"
        assert annotation_version.genome_build == vcf.genome_build, msg
        calculate_cohort_stats(vcf.cohort, annotation_version)
        _compute_vcf_specific_stats(vcf, annotation_version)
    except:
        tb = get_traceback()
        create_event(None, 'calculate_vcf_stats', details=tb, severity=LogLevel.ERROR)
        raise


def _get_sample_stats_code_version() -> SampleStatsCodeVersion:
    """
        Version | Notes on sample stats code changes
        0 - Set for historical versions created before code_version added
        1 - Existing code when code_version added
        2 - Added VCFLengthStats
        3 - Optimisations
        4 - fusions_count, and gene-level variants left out of VCFLengthStats
        5 - Per-sample totals and chrX hom only count the sample's own variant calls (#2087)
    """
    code_version, _ = thread_safe_unique_together_get_or_create(SampleStatsCodeVersion,
                                                                name="SampleStats",
                                                                version=SAMPLE_STATS_CODE_VERSION,
                                                                code_version=CodeVersion.current())
    return code_version


def _cohort_stats_are_fresh(cgc, annotation_version, code_version) -> bool:
    """ True when this CGC already has a genotype-level per-sample row at the current
        code_version plus matching variant/gene/clinvar annotation-version rows, i.e.
        a recompute would reproduce what is already stored. """
    vav = annotation_version.variant_annotation_version
    gav = annotation_version.gene_annotation_version
    cv = annotation_version.clinvar_version
    return (
        CohortGenotypeStats.objects.filter(
            cohort_genotype_collection=cgc, code_version=code_version).exists()
        and CohortGenotypeVariantAnnotationStats.objects.filter(
            cohort_genotype_collection=cgc, code_version=code_version,
            variant_annotation_version=vav).exists()
        and CohortGenotypeGeneAnnotationStats.objects.filter(
            cohort_genotype_collection=cgc, code_version=code_version,
            gene_annotation_version=gav).exists()
        and CohortGenotypeClinVarAnnotationStats.objects.filter(
            cohort_genotype_collection=cgc, code_version=code_version,
            clinvar_version=cv).exists()
    )


def _compute_vcf_specific_stats(vcf: VCF, annotation_version: AnnotationVersion):
    """ VCF-keyed stats — VEP-skipped count and variant-class length histograms.
        Per-sample / cohort-aggregate stats are written by calculate_cohort_stats. """
    start = time.time()
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs = qs.filter(Variant.get_no_reference_q())
    qs = qs.exclude(Variant.get_gene_level_q())  # no coordinate, so no length
    qs = vcf.get_variant_qs(qs)
    vep_skipped_count = qs.filter(variantannotation__vep_skipped_reason__isnull=False) \
        .exclude(variantannotation__vep_skipped_reason="").count()

    # SNVs have no length histogram and are most of the rows, so leave them in the database
    non_snv_qs = qs.annotate(ref_len=Length("locus__ref__seq"), alt_len=Length("alt__seq")) \
        .exclude(ref_len=1, alt_len=1)
    columns = ["alt__seq", "svlen", "variantannotation__hgvs_g", "ref_len", "alt_len"]
    variant_class_lengths = defaultdict(list)
    for alt, svlen, hgvs_g, ref_len, alt_len in non_snv_qs.values_list(*columns).iterator(chunk_size=10_000):
        variant_class = _get_variant_class(alt, hgvs_g)
        variant_length = abs(svlen) if svlen else max(alt_len, ref_len)
        variant_class_lengths[variant_class].append(variant_length)

    if vep_skipped_count:
        stats, created = VCFAnnotationStats.objects.get_or_create(
            vcf=vcf, variant_annotation_version=annotation_version.variant_annotation_version,
            defaults={"vep_skipped_count": vep_skipped_count})
        if not created and vep_skipped_count > stats.vep_skipped_count:
            stats.vep_skipped_count = vep_skipped_count
            stats.save()

    if variant_class_lengths:
        code_version = _get_sample_stats_code_version()
        collection, created = VCFLengthStatsCollection.objects.update_or_create(
            vcf=vcf, defaults={"code_version": code_version})
        if not created:
            collection.vcflengthstats_set.all().delete()

        vcf_length_stats = []
        for variant_class, lengths in variant_class_lengths.items():
            lengths = np.array(lengths)
            is_log = False
            abs_lengths = np.abs(lengths)
            max_min_ratio = abs_lengths.max() / (abs_lengths.min() + 1)
            logging.info(f"Variant class {variant_class} has {abs_lengths.max()=} / {abs_lengths.min()=} = {max_min_ratio=}")
            if max_min_ratio > 1000:  # Spans 3 orders of magnitude
                is_log = True
                lengths = np.log10(lengths)

            histogram_counts, histogram_bin_edges = np.histogram(lengths, bins=50)
            vcf_length_stats.append(VCFLengthStats(
                collection=collection, variant_class=variant_class, is_log=is_log,
                histogram_counts=histogram_counts.tolist(),
                histogram_bin_edges=histogram_bin_edges.tolist()))
        if vcf_length_stats:
            VCFLengthStats.objects.bulk_create(vcf_length_stats)

    end = time.time()
    logging.info("VCF-specific stats for %s took %.2f seconds", vcf, end - start)


@celery.shared_task
def calculate_cohort_stats_task(cohort_id, annotation_version_id):
    cohort = Cohort.objects.get(pk=cohort_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    calculate_cohort_stats(cohort, annotation_version)


# Field indexes for the per-bucket counter arrays. Keep in lock-step with
# _COUNTER_FIELDS below — index N in the array is the field named at position N.
_COUNTER_FIELDS = (
    "variant_count",
    "snp_count",
    "insertions_count",
    "deletions_count",
    "fusions_count",
    "ref_count",
    "het_count",
    "hom_count",
    "unk_count",
    "x_hom_count",
    "x_het_count",
    "x_unk_count",
    "variant_dbsnp_count",
    "snp_dbsnp_count",
    "insertions_dbsnp_count",
    "deletions_dbsnp_count",
    "ref_high_or_moderate_count",
    "het_high_or_moderate_count",
    "hom_high_or_moderate_count",
    "unk_high_or_moderate_count",
    "gene_count",
    "ref_omim_phenotype_count",
    "het_omim_phenotype_count",
    "hom_omim_phenotype_count",
    "unk_omim_phenotype_count",
    "clinvar_count",
    "ref_clinvar_pathogenic_count",
    "het_clinvar_pathogenic_count",
    "hom_clinvar_pathogenic_count",
    "unk_clinvar_pathogenic_count",
)
_FIELD_INDEX = {name: i for i, name in enumerate(_COUNTER_FIELDS)}
_NUM_FIELDS = len(_COUNTER_FIELDS)


def _new_counter():
    return [0] * _NUM_FIELDS


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


def get_filter_keys_to_precompute_for_cohort(cohort) -> list[Optional[str]]:
    """ The filter_key buckets calculate_cohort_stats populates: the raw aggregate (None) row, plus
        the inheritance modes when the cohort has a Trio. """
    keys = [None]
    if trio := cohort.trio_set.first():
        keys.extend(trio.get_stats_filter_keys())
    return keys


def calculate_cohort_stats(cohort: Cohort, annotation_version: AnnotationVersion):
    """ Count the cohort's rows once (@see _count_cohort_stats) and populate the four
        CohortGenotype*Stats classes for:
          - per-sample rows (sample IS NOT NULL, filter_key NULL) — every cohort member
          - aggregate row (sample IS NULL, filter_key NULL)
          - filter-keyed aggregate rows (sample IS NULL, filter_key=...) for any keys
            returned by get_filter_keys_to_precompute_for_cohort(cohort)
        All rows are FK'd to cohort.cohort_genotype_collection (the current CGC).
        Idempotent: any existing rows for this CGC + the four annotation versions
        are deleted and rewritten.

        Serialized per CGC: the CohortGenotypeCollection row is locked with
        select_for_update() inside a transaction so two concurrent recomputes for
        the same CGC can't race the delete/bulk_create (issue #1576). A task that
        acquires the lock after another has just finished early-exits via
        _cohort_stats_are_fresh instead of redoing the work. Unrelated cohorts
        lock different rows and still run concurrently. """
    code_version = _get_sample_stats_code_version()
    with transaction.atomic():
        cgc = CohortGenotypeCollection.objects.select_for_update().get(
            pk=cohort.cohort_genotype_collection.pk)
        if _cohort_stats_are_fresh(cgc, annotation_version, code_version):
            logging.info("Cohort stats already fresh for %s; skipping recompute", cgc)
            return
        _compute_and_persist_cohort_stats(cohort, cgc, annotation_version, code_version)


def _compute_and_persist_cohort_stats(cohort: Cohort, cgc: CohortGenotypeCollection,
                                      annotation_version: AnnotationVersion,
                                      code_version: SampleStatsCodeVersion):
    """ Count the four CohortGenotype*Stats families for the cohort and persist
        them via _persist_cohort_stats. Must run inside the per-CGC locked
        transaction opened by calculate_cohort_stats, using the locked cgc. """
    has_filters = _cohort_has_filters(cohort)
    precompute_keys = get_filter_keys_to_precompute_for_cohort(cohort)
    sample_counters, aggregate_counters = _count_cohort_stats(cohort, cgc, annotation_version,
                                                              has_filters, precompute_keys)
    samples = [cs.sample for cs in cohort.get_cohort_samples()]
    _persist_cohort_stats(cgc, annotation_version, code_version, has_filters,
                          sample_counters, aggregate_counters, precompute_keys, samples)


class _VariantFlags(NamedTuple):
    """ What the counters need to know about a cohort row, besides its zygosities. Rows sharing a
        samples_zygosity string and these flags are counted once in SQL and folded in together """
    symbolic_alt: Optional[str]
    ref_alt_length_sign: int  # sign(len(ref) - len(alt))
    is_snp: bool
    dbsnp: bool
    high_or_moderate: bool
    has_transcript: bool
    omim: bool
    clinvar: bool
    clinvar_pathogenic: bool
    is_x: bool
    passing_filter: bool


# SQL for each boolean _VariantFlags field over the columns of _get_cohort_rows_queryset, written to match
# Python truthiness: an empty dbsnp/transcript/OMIM string or a highest_pathogenicity of 0 is False
_BOOLEAN_FLAGS_SQL = {
    "is_snp": "length(v.ref_seq) = 1 AND length(v.alt_seq) = 1",
    "dbsnp": "coalesce(v.dbsnp_rs_id, '') <> ''",
    "high_or_moderate": "coalesce(v.impact IN (%s, %s), false)",
    "has_transcript": "coalesce(v.transcript_id, '') <> ''",
    "omim": "coalesce(v.omim_terms, '') <> ''",
    "clinvar": "coalesce(v.clinvar_highest_pathogenicity, 0) <> 0",
    "clinvar_pathogenic": "coalesce(v.clinvar_highest_pathogenicity, 0) >= %s",
    "is_x": "v.contig = %s",
    "passing_filter": "coalesce(v.cohort_filters, '') = ''",
}
_BOOLEAN_FLAGS_PARAMS = [PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE, 4, "X"]
assert tuple(_BOOLEAN_FLAGS_SQL) == _VariantFlags._fields[2:]

# Grouping on one integer is much cheaper than on each flag. Bit i is boolean flag i, the bits above
# them are ref_alt_length_sign + 1
_LENGTH_SIGN_SHIFT = len(_BOOLEAN_FLAGS_SQL)
_FLAG_BITS_SQL = " | ".join([f"CASE WHEN {sql} THEN {1 << i} ELSE 0 END" for i, sql in enumerate(_BOOLEAN_FLAGS_SQL.values())]
                            + [f"((sign(length(v.ref_seq) - length(v.alt_seq))::integer + 1) << {_LENGTH_SIGN_SHIFT})"])
_SYMBOLIC_ALT_SQL = "CASE WHEN v.alt_seq LIKE '<%%>' THEN v.alt_seq END"


@cache
def _decode_flags(symbolic_alt: Optional[str], flag_bits: int) -> _VariantFlags:
    booleans = [bool(flag_bits & (1 << i)) for i in range(_LENGTH_SIGN_SHIFT)]
    return _VariantFlags(symbolic_alt, (flag_bits >> _LENGTH_SIGN_SHIFT) - 1, *booleans)


# Fields a zygosity counts into: total, high/moderate impact, OMIM, ClinVar pathogenic, chrX.
# chrX hom/het is the sex guess (@see CohortGenotypeStats.chrx_sex_guess). A hom ref call says nothing
# about hemizygosity, and in a joint-called VCF it's another sample's variant, so it has no chrX field
_ZYGOSITY_FIELDS = {
    Zygosity.HOM_ALT: ("hom_count", "hom_high_or_moderate_count", "hom_omim_phenotype_count",
                       "hom_clinvar_pathogenic_count", "x_hom_count"),
    Zygosity.HET: ("het_count", "het_high_or_moderate_count", "het_omim_phenotype_count",
                   "het_clinvar_pathogenic_count", "x_het_count"),
    Zygosity.HOM_REF: ("ref_count", "ref_high_or_moderate_count", "ref_omim_phenotype_count",
                       "ref_clinvar_pathogenic_count", None),
    Zygosity.UNKNOWN_ZYGOSITY: ("unk_count", "unk_high_or_moderate_count", "unk_omim_phenotype_count",
                                "unk_clinvar_pathogenic_count", "x_unk_count"),
}
_ZYGOSITY_FIELD_INDEXES = {
    zygosity: tuple(None if f is None else _FIELD_INDEX[f] for f in fields)
    for zygosity, fields in _ZYGOSITY_FIELDS.items()
}

# The aggregate classifies a row by the first of these in its zygosity string, else unknown
_AGGREGATE_ZYGOSITY_PRIORITY = (Zygosity.HOM_ALT, Zygosity.HET, Zygosity.HOM_REF)

_SYMBOLIC_INSERTIONS = {VCFSymbolicAllele.DUP, VCFSymbolicAllele.INS}
_FUSION_KINDS = {GeneLevelSymbolicAlt.FUSION, GeneLevelSymbolicAlt.FUSION_UNORDERED}


def _variant_class_field(flags: _VariantFlags) -> str:
    if flags.is_snp:
        return "snp_count"
    if alt := flags.symbolic_alt:
        if alt in _SYMBOLIC_INSERTIONS:
            return "insertions_count"
        if alt == VCFSymbolicAllele.DEL:
            return "deletions_count"
        if (gene_level := GeneLevelSymbolicAlt.parse(alt)) and gene_level[0] in _FUSION_KINDS:
            return "fusions_count"
        return "snp_count"  # fallback
    if flags.ref_alt_length_sign > 0:
        return "deletions_count"
    return "insertions_count"


def _dbsnp_class_field(flags: _VariantFlags) -> Optional[str]:
    if flags.is_snp:
        return "snp_dbsnp_count"
    if flags.ref_alt_length_sign > 0:
        return "deletions_dbsnp_count"
    if flags.ref_alt_length_sign < 0:
        return "insertions_dbsnp_count"
    return None


@cache
def _flag_field_indexes(flags: _VariantFlags) -> tuple[list[int], list[int]]:
    """ (fields a row with a call adds to, impact/OMIM/ClinVar slots of _ZYGOSITY_FIELDS it adds to) """
    totals = [_FIELD_INDEX["variant_count"], _FIELD_INDEX[_variant_class_field(flags)]]
    if flags.dbsnp:
        totals.append(_FIELD_INDEX["variant_dbsnp_count"])
        if dbsnp_class_field := _dbsnp_class_field(flags):
            totals.append(_FIELD_INDEX[dbsnp_class_field])
    if flags.has_transcript:
        totals.append(_FIELD_INDEX["gene_count"])
    if flags.clinvar:
        totals.append(_FIELD_INDEX["clinvar_count"])
    zygosity_slots = [0]
    if flags.high_or_moderate:
        zygosity_slots.append(1)
    if flags.omim:
        zygosity_slots.append(2)
    if flags.clinvar_pathogenic:
        zygosity_slots.append(3)
    return totals, zygosity_slots


def _add_count(counter: list[int], flags: _VariantFlags, zygosity: str, n: int,
               has_call: bool = True, count_chrx: bool = True):
    """ Add n rows with these flags, classified under zygosity, to a counter array. Rows without a
        call still count in the per-zygosity tallies - SampleNode reads them as cached counts for its
        zygosity checkboxes (@see stats_cache.get_cached_label_count_for_cohort) """
    totals, zygosity_slots = _flag_field_indexes(flags)
    if has_call:
        for i in totals:
            counter[i] += n
    zygosity_fields = _ZYGOSITY_FIELD_INDEXES.get(zygosity, _ZYGOSITY_FIELD_INDEXES[Zygosity.UNKNOWN_ZYGOSITY])
    for slot in zygosity_slots:
        counter[zygosity_fields[slot]] += n
    if count_chrx and flags.is_x and (x_field := zygosity_fields[4]) is not None:
        counter[x_field] += n


def _get_cohort_rows_queryset(cgc: CohortGenotypeCollection, annotation_version: AnnotationVersion,
                              with_gene: bool):
    """ One row per cohort variant with the columns _BOOLEAN_FLAGS_SQL reads, through the annotation
        version's partitions and every collection the CGC reads (a sub-cohort reads its parent's) """
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs = qs.filter(Variant.get_no_reference_q())
    qs = qs.annotate(**cgc.get_annotation_kwargs())
    qs = qs.filter(**{f"{cgc.cohortgenotype_alias}__isnull": False})
    columns = {
        "contig": "locus__contig__name",
        "ref_seq": "locus__ref__seq",
        "alt_seq": "alt__seq",
        "samples_zygosity": f"{cgc.cohortgenotype_alias}__samples_zygosity",
        "dbsnp_rs_id": "variantannotation__dbsnp_rs_id",
        "impact": "variantannotation__impact",
        "transcript_id": "variantannotation__transcript_id",
        "omim_terms": "variantannotation__gene__geneannotation__omim_terms",
        "clinvar_highest_pathogenicity": "clinvar__highest_pathogenicity",
        "cohort_filters": f"{cgc.cohortgenotype_alias}__filters",
    }
    if with_gene:
        columns["gene"] = "variantannotation__gene_id"
    return qs.values(**{alias: F(path) for alias, path in columns.items()})


def _count_cohort_stats(cohort: Cohort, cgc: CohortGenotypeCollection, annotation_version: AnnotationVersion,
                        has_filters: bool, precompute_keys: list[Optional[str]]):
    """ Returns (sample_counters, aggregate_counters) for _persist_cohort_stats:
            sample_counters: (sample_pk, passing_filter) → counter array
            aggregate_counters: (filter_key, passing_filter) → counter array

        Postgres counts: the cohort rows are grouped by (samples_zygosity, _VariantFlags) into a temp
        table - many rows share a zygosity string - then grouped three ways: per packed sample index,
        by the aggregate's priority zygosity, and by the trio's three zygosities. Python folds the
        small results onto the counter arrays, and keeps the trio predicates (shared with
        trio_node.py) and the compound het gene pairing (#2283).
        Must run inside a transaction - the temp table is ON COMMIT DROP """
    cohort_samples = list(cohort.get_cohort_samples())
    trio, trio_indexes = _get_trio_and_packed_indexes(cohort, cohort_samples)
    sample_counters: dict[tuple[int, bool], list[int]] = defaultdict(_new_counter)
    aggregate_counters: dict[tuple[Optional[str], bool], list[int]] = defaultdict(_new_counter)

    def passing_filter_buckets(flags: _VariantFlags):
        return (False, True) if has_filters and flags.passing_filter else (False,)

    with_gene = trio_indexes is not None
    rows_sql, rows_params = _get_cohort_rows_queryset(cgc, annotation_version, with_gene).query.sql_with_params()
    table = f"pg_temp.cohort_stats_cgc_{cgc.pk}"
    gene_column = ", v.gene" if with_gene else ""
    group_by_positions = "1, 2, 3, 4" if with_gene else "1, 2, 3"

    # The temp table's GROUP BY spills to disk at the default work_mem on wide cohorts, where almost every
    # zygosity string is distinct
    with (pg_settings(local=True, jit="off", work_mem=settings.ANNOTATION_COHORT_STATS_WORK_MEM),
          connection.cursor() as cursor):
        cursor.execute(f"CREATE TEMPORARY TABLE {table} ON COMMIT DROP AS "
                       f"SELECT v.samples_zygosity, {_SYMBOLIC_ALT_SQL} AS symbolic_alt, "
                       f"{_FLAG_BITS_SQL} AS flag_bits{gene_column}, count(*) AS n "
                       f"FROM ({rows_sql}) v GROUP BY {group_by_positions}",
                       [*_BOOLEAN_FLAGS_PARAMS, *rows_params])

        # Per-sample: string_to_table(.., NULL) splits into characters - substr per index would walk the
        # string from the start each time, O(samples²) per row. Rows split by gene for the trio are merged
        # first, as the split almost doubles the rows of a wide cohort. A sub-cohort reads its parent's
        # strings, so the other samples' characters are dropped before they are grouped
        without_gene = f"(SELECT samples_zygosity, symbolic_alt, flag_bits, sum(n) AS n FROM {table} GROUP BY 1, 2, 3)" \
            if with_gene else table
        sub_cohort_where = ""
        sub_cohort_params = []
        if cohort.is_sub_cohort:
            sub_cohort_where = "WHERE s.packed_index = ANY(%s)"
            sub_cohort_params = [[cs.cohort_genotype_packed_field_index + 1 for cs in cohort_samples]]
        cursor.execute(f"SELECT s.packed_index - 1, s.zygosity, g.symbolic_alt, g.flag_bits, sum(g.n)::bigint "
                       f"FROM {without_gene} g CROSS JOIN LATERAL "
                       f"string_to_table(g.samples_zygosity, NULL) WITH ORDINALITY AS s(zygosity, packed_index) "
                       f"{sub_cohort_where} GROUP BY 1, 2, 3, 4", sub_cohort_params)
        # A sample's totals count the rows it has a variant call on - in a multi-sample VCF the other rows
        # are the other samples' variants. Without a GT field every call is unknown zygosity, so count them all.
        sample_by_packed_index = {}
        for cs in cohort_samples:
            call_zygosities = {Zygosity.HET, Zygosity.HOM_ALT}
            if not cs.sample.has_genotype:
                call_zygosities.add(Zygosity.UNKNOWN_ZYGOSITY)
            sample_by_packed_index[cs.cohort_genotype_packed_field_index] = (cs.sample.pk, call_zygosities)
        for packed_index, zygosity, symbolic_alt, flag_bits, n in cursor.fetchall():
            if sample_and_calls := sample_by_packed_index.get(packed_index):
                sample_pk, call_zygosities = sample_and_calls
                flags = _decode_flags(symbolic_alt, flag_bits)
                for pf_pass in passing_filter_buckets(flags):
                    _add_count(sample_counters[(sample_pk, pf_pass)], flags, zygosity, n,
                               has_call=zygosity in call_zygosities)

        # Aggregate: every cohort row counts, as a CohortNode with no zygosity filter returns every row,
        # classified by the highest-priority zygosity present across the cohort's samples - hom > het > ref > unk.
        # This makes ref+het+hom+unk_count sum to variant_count, the invariant CohortGenotypeStats relies on.
        priority_sql = " ".join("WHEN strpos(g.samples_zygosity, %s) > 0 THEN %s" for _ in _AGGREGATE_ZYGOSITY_PRIORITY)
        priority_params = [z for zygosity in _AGGREGATE_ZYGOSITY_PRIORITY for z in (zygosity, zygosity)]
        cursor.execute(f"SELECT CASE {priority_sql} ELSE %s END, g.symbolic_alt, g.flag_bits, sum(g.n)::bigint "
                       f"FROM {table} g GROUP BY 1, 2, 3",
                       [*priority_params, Zygosity.UNKNOWN_ZYGOSITY])
        for zygosity, symbolic_alt, flag_bits, n in cursor.fetchall():
            flags = _decode_flags(symbolic_alt, flag_bits)
            for pf_pass in passing_filter_buckets(flags):
                _add_count(aggregate_counters[(None, pf_pass)], flags, zygosity, n)

        if with_gene:
            # Only a proband HET can be compound het (_compound_het_source), so only those rows keep their gene
            mum, dad, proband = (f"substr(g.samples_zygosity, {i + 1}, 1)" for i in trio_indexes)  # 1-based
            cursor.execute(f"SELECT {mum}, {dad}, {proband}, CASE WHEN {proband} = %s THEN g.gene END, "
                           f"g.symbolic_alt, g.flag_bits, sum(g.n)::bigint "
                           f"FROM {table} g GROUP BY 1, 2, 3, 4, 5, 6",
                           [Zygosity.HET])
            _count_trio_modes(trio, cursor.fetchall(), precompute_keys, aggregate_counters, passing_filter_buckets)
        cursor.execute(f"DROP TABLE {table}")
    return sample_counters, aggregate_counters


def _count_trio_modes(trio: Trio, trio_rows, precompute_keys, aggregate_counters, passing_filter_buckets):
    """ trio_rows are (mum, dad, proband zygosity, gene_id, symbolic_alt, flag_bits, n). Each inheritance mode
        bucket counts like the aggregate, classified by the proband's zygosity, without chrX fields.
        Compound het rows (proband HET from exactly one parent) only count for genes that have a
        mum-only and a dad-only source """
    chet_key = trio.get_stats_filter_key("compound_het")
    chet_gene_sources: dict[tuple[int, bool], set] = defaultdict(set)
    chet_pending = []
    for mum_zyg, dad_zyg, proband_zyg, gene_id, symbolic_alt, flag_bits, n in trio_rows:
        flags = _decode_flags(symbolic_alt, flag_bits)
        modes = _trio_predicates(mum_zyg, dad_zyg, proband_zyg, flags.is_x,
                                 trio.mother_affected, trio.father_affected)
        chet_source = _compound_het_source(mum_zyg, dad_zyg, proband_zyg)
        for pf_pass in passing_filter_buckets(flags):
            for mode in modes:
                key = trio.get_stats_filter_key(mode)
                if key in precompute_keys:
                    _add_count(aggregate_counters[(key, pf_pass)], flags, proband_zyg, n, count_chrx=False)
            if chet_source is not None and gene_id is not None:
                chet_gene_sources[(gene_id, pf_pass)].add(chet_source)
                chet_pending.append((gene_id, pf_pass, flags, proband_zyg, n))

    for gene_id, pf_pass, flags, proband_zyg, n in chet_pending:
        if {"mum", "dad"} <= chet_gene_sources[(gene_id, pf_pass)]:
            _add_count(aggregate_counters[(chet_key, pf_pass)], flags, proband_zyg, n, count_chrx=False)


def _compound_het_source(mum_zyg, dad_zyg, proband_zyg) -> Optional[str]:
    """ The parent a proband HET is inherited from, when exactly one parent is HET """
    no_variant = (Zygosity.MISSING, Zygosity.HOM_REF)
    if proband_zyg == Zygosity.HET:
        if mum_zyg == Zygosity.HET and dad_zyg in no_variant:
            return "mum"
        if mum_zyg in no_variant and dad_zyg == Zygosity.HET:
            return "dad"
    return None


def _get_trio_and_packed_indexes(cohort: Cohort, cohort_samples) -> tuple[Optional[Trio], Optional[tuple[int, int, int]]]:
    """ The cohort's trio and its (mother, father, proband) packed indexes, or None indexes when the
        trio isn't coherent with this cohort, which skips inheritance bucketing """
    trio = cohort.trio_set.first()
    if trio is None:
        return None, None
    packed_index_by_cohort_sample_id = {cs.pk: cs.cohort_genotype_packed_field_index for cs in cohort_samples}
    try:
        return trio, (packed_index_by_cohort_sample_id[trio.mother_id],
                      packed_index_by_cohort_sample_id[trio.father_id],
                      packed_index_by_cohort_sample_id[trio.proband_id])
    except KeyError:
        return trio, None


def _cohort_has_filters(cohort: Cohort) -> bool:
    """ A cohort has 'filters' iff its underlying VCF (if any) does. Custom
        multi-VCF cohorts inherit filtering from their member VCFs but, for
        passing-filter stats purposes, we mirror the per-VCF flag of the
        cohort's base VCF if any. """
    vcf = cohort.get_vcf()
    return bool(vcf and vcf.has_filters)


def _persist_cohort_stats(cgc, annotation_version, code_version, has_filters,
                          sample_counters, aggregate_counters, precompute_keys, samples):
    """ Replace any existing rows for this CGC + annotation versions on the
        new models with the computed counters. """
    vav = annotation_version.variant_annotation_version
    gav = annotation_version.gene_annotation_version
    cv = annotation_version.clinvar_version

    # Wipe existing rows for this CGC scoped by annotation version where applicable
    CohortGenotypeStats.objects.filter(cohort_genotype_collection=cgc).delete()
    CohortGenotypeVariantAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, variant_annotation_version=vav).delete()
    CohortGenotypeGeneAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, gene_annotation_version=gav).delete()
    CohortGenotypeClinVarAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, clinvar_version=cv).delete()

    rows_genotype = []
    rows_variant = []
    rows_gene = []
    rows_clinvar = []

    def _populate(target, counter, fields):
        for f in fields:
            setattr(target, f, counter[_FIELD_INDEX[f]])

    GENOTYPE_FIELDS = (
        "variant_count", "snp_count", "insertions_count", "deletions_count", "fusions_count",
        "ref_count", "het_count", "hom_count", "unk_count",
        "x_hom_count", "x_het_count", "x_unk_count",
    )
    VARIANT_ANN_FIELDS = (
        "variant_dbsnp_count", "snp_dbsnp_count", "insertions_dbsnp_count", "deletions_dbsnp_count",
        "ref_high_or_moderate_count", "het_high_or_moderate_count",
        "hom_high_or_moderate_count", "unk_high_or_moderate_count",
    )
    GENE_ANN_FIELDS = (
        "gene_count",
        "ref_omim_phenotype_count", "het_omim_phenotype_count",
        "hom_omim_phenotype_count", "unk_omim_phenotype_count",
    )
    CLINVAR_ANN_FIELDS = (
        "clinvar_count",
        "ref_clinvar_pathogenic_count", "het_clinvar_pathogenic_count",
        "hom_clinvar_pathogenic_count", "unk_clinvar_pathogenic_count",
    )

    sample_by_pk = {s.pk: s for s in samples}

    # Per-sample rows (filter_key NULL)
    for (sample_pk, pf_pass), counter in sample_counters.items():
        if pf_pass and not has_filters:
            continue
        sample = sample_by_pk[sample_pk]
        gs = CohortGenotypeStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
        )
        _populate(gs, counter, GENOTYPE_FIELDS)
        rows_genotype.append(gs)
        va = CohortGenotypeVariantAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
        )
        _populate(va, counter, VARIANT_ANN_FIELDS)
        rows_variant.append(va)
        ga = CohortGenotypeGeneAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
        )
        _populate(ga, counter, GENE_ANN_FIELDS)
        rows_gene.append(ga)
        ca = CohortGenotypeClinVarAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
        )
        _populate(ca, counter, CLINVAR_ANN_FIELDS)
        rows_clinvar.append(ca)

    # Aggregate / filter-keyed rows (sample NULL)
    for (filter_key, pf_pass), counter in aggregate_counters.items():
        if pf_pass and not has_filters:
            continue
        gs = CohortGenotypeStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
        )
        _populate(gs, counter, GENOTYPE_FIELDS)
        rows_genotype.append(gs)
        va = CohortGenotypeVariantAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
        )
        _populate(va, counter, VARIANT_ANN_FIELDS)
        rows_variant.append(va)
        ga = CohortGenotypeGeneAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
        )
        _populate(ga, counter, GENE_ANN_FIELDS)
        rows_gene.append(ga)
        ca = CohortGenotypeClinVarAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
        )
        _populate(ca, counter, CLINVAR_ANN_FIELDS)
        rows_clinvar.append(ca)

    passing_filter_values = (False, True) if has_filters else (False,)

    # Ensure every sample has a row even if it had no variants, otherwise readers
    # treat the absence as stale and re-enqueue a recompute on every view.
    for sample in samples:
        for pf_pass in passing_filter_values:
            if (sample.pk, pf_pass) in sample_counters:
                continue
            rows_genotype.append(CohortGenotypeStats(
                cohort_genotype_collection=cgc, sample=sample, filter_key=None,
                passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
            ))
            rows_variant.append(CohortGenotypeVariantAnnotationStats(
                cohort_genotype_collection=cgc, sample=sample, filter_key=None,
                passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
            ))
            rows_gene.append(CohortGenotypeGeneAnnotationStats(
                cohort_genotype_collection=cgc, sample=sample, filter_key=None,
                passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
            ))
            rows_clinvar.append(CohortGenotypeClinVarAnnotationStats(
                cohort_genotype_collection=cgc, sample=sample, filter_key=None,
                passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
            ))

    # Ensure every precomputed (key, pf_pass) bucket has a row even if no
    # variants matched it (so reads hit instead of falling back to live count).
    existing_agg_keys = {(fk, pf) for fk, pf in aggregate_counters.keys()}
    for fk in precompute_keys:
        for pf_pass in passing_filter_values:
            if (fk, pf_pass) in existing_agg_keys:
                continue
            rows_genotype.append(CohortGenotypeStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
            ))
            rows_variant.append(CohortGenotypeVariantAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
            ))
            rows_gene.append(CohortGenotypeGeneAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
            ))
            rows_clinvar.append(CohortGenotypeClinVarAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
            ))

    CohortGenotypeStats.objects.bulk_create(rows_genotype, batch_size=2000)
    CohortGenotypeVariantAnnotationStats.objects.bulk_create(rows_variant, batch_size=2000)
    CohortGenotypeGeneAnnotationStats.objects.bulk_create(rows_gene, batch_size=2000)
    CohortGenotypeClinVarAnnotationStats.objects.bulk_create(rows_clinvar, batch_size=2000)


def _get_variant_class(alt: str, hgvs_g: str) -> Optional[VariantClass]:
    """ The VEP VariantClass is not always accurate so basically DIY here """

    variant_type = None
    if Sequence.allele_is_symbolic(alt):
        if alt == VCFSymbolicAllele.INV:
            variant_type = VariantClass.INVERSION
        elif alt == VCFSymbolicAllele.DEL:
            variant_type = VariantClass.DELETION
        elif alt == VCFSymbolicAllele.DUP:
            variant_type = VariantClass.DUPLICATION
    elif hgvs_g:
        if "delins" in hgvs_g:
            variant_type = VariantClass.INDEL
        elif "inv" in hgvs_g:
            variant_type = VariantClass.INVERSION
        elif "ins" in hgvs_g:
            variant_type = VariantClass.INSERTION
        elif "del" in hgvs_g:
            variant_type = VariantClass.DELETION
        elif "dup" in hgvs_g:
            variant_type = VariantClass.DUPLICATION
    return variant_type


def enqueue_cohort_stats_recompute(cohort: Cohort, annotation_version: AnnotationVersion) -> None:
    """ Fire-and-forget recompute for readers that missed the cache. Idempotent —
        calculate_cohort_stats early-exits if the rows are already fresh. """
    calculate_cohort_stats_task.apply_async(
        args=[cohort.pk, annotation_version.pk],
        queue="annotation_workers",
    )
