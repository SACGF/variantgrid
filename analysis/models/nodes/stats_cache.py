"""
    Reader-side cache lookups for the CohortGenotype*Stats family. Each source node translates its
    config to a filter_key in get_stats_filter_key.

    The writer side (annotation/tasks/calculate_sample_stats.calculate_cohort_stats) and the readers
    both build inheritance keys through Trio.get_stats_filter_key. A silent mismatch between the two
    would produce cache misses forever.
"""
from typing import Optional, Union

from django.core.exceptions import ObjectDoesNotExist

from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from annotation.tasks.calculate_sample_stats import enqueue_cohort_stats_recompute
from snpdb.models import CohortGenotypeStats
from snpdb.models.models_enums import BuiltInFilters


class _Uncacheable:
    """ Sentinel returned by get_stats_filter_key for node configs we don't precompute. """


UNCACHEABLE = _Uncacheable()
FilterKey = Union[None, str, _Uncacheable]


LABEL_TO_MODEL = {
    BuiltInFilters.TOTAL: CohortGenotypeStats,
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: CohortGenotypeVariantAnnotationStats,
    BuiltInFilters.OMIM: CohortGenotypeGeneAnnotationStats,
    BuiltInFilters.CLINVAR: CohortGenotypeClinVarAnnotationStats,
}


def get_cached_label_count_for_cohort(cohort, sample, filter_key: FilterKey,
                                      annotation_version, passing_filter,
                                      zygosities, label) -> Optional[int]:
    """ sample=None → aggregate row; sample=Sample(...) → per-sample row.
        filter_key=None → raw count; canonical-JSON string → filter-keyed row.
        Returns None on cache miss (caller falls back to live count); on miss
        also fires a lazy recompute task. """
    if filter_key is UNCACHEABLE:
        return None
    model = LABEL_TO_MODEL.get(label)
    if model is None:
        return None
    cgc = cohort.cohort_genotype_collection
    lookup = {
        "cohort_genotype_collection": cgc,
        "sample": sample,
        "filter_key": filter_key,
        "passing_filter": passing_filter,
    }
    if model is CohortGenotypeStats:
        pass
    elif model is CohortGenotypeVariantAnnotationStats:
        lookup["variant_annotation_version"] = annotation_version.variant_annotation_version
    elif model is CohortGenotypeGeneAnnotationStats:
        lookup["gene_annotation_version"] = annotation_version.gene_annotation_version
    elif model is CohortGenotypeClinVarAnnotationStats:
        lookup["clinvar_version"] = annotation_version.clinvar_version

    try:
        obj = model.objects.get(**lookup)
    except ObjectDoesNotExist:
        _enqueue_lazy_recompute(cohort, annotation_version)
        return None
    return obj.count_for_zygosity(*zygosities, label=label)


def _enqueue_lazy_recompute(cohort, annotation_version) -> None:
    """ Fire-and-forget. Idempotent — task no-ops if the row now exists. """
    enqueue_cohort_stats_recompute(cohort, annotation_version)
