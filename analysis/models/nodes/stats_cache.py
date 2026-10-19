"""
    Reader-side cache lookups for the CohortGenotype*Stats family. A source node that can read the
    cache passes the filter_key for its config: None for the raw aggregate row, or an inheritance
    bucket.

    The writer side (annotation/tasks/calculate_sample_stats.calculate_cohort_stats) and the readers
    both build keys through Cohort.get_stats_filter_key and Trio.get_stats_filter_key. A silent mismatch
    between the two would produce cache misses forever. A sub cohort reads its parent's CGC: its
    per-sample rows are the parent's, its aggregate rows are keyed by its samples.
"""
from typing import Optional

from django.core.exceptions import ObjectDoesNotExist

from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from annotation.tasks.calculate_sample_stats import enqueue_cohort_stats_recompute
from snpdb.models import CohortGenotypeStats
from snpdb.models.models_enums import BuiltInFilters

LABEL_TO_MODEL = {
    BuiltInFilters.TOTAL: CohortGenotypeStats,
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: CohortGenotypeVariantAnnotationStats,
    BuiltInFilters.OMIM: CohortGenotypeGeneAnnotationStats,
    BuiltInFilters.CLINVAR: CohortGenotypeClinVarAnnotationStats,
}


def get_cached_label_count_for_cohort(cohort, sample, filter_key: Optional[str],
                                      annotation_version, passing_filter,
                                      zygosities, label) -> Optional[int]:
    """ sample=None → aggregate row; sample=Sample(...) → per-sample row.
        filter_key=None → raw count; canonical-JSON string → filter-keyed row.
        Returns None on cache miss (caller falls back to live count); on miss
        also fires a lazy recompute task. """
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
