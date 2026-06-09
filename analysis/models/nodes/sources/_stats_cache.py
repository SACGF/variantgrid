"""
    Reader-side cache lookups for the CohortGenotype*Stats family + node-class
    handlers that translate node configs to canonical filter_keys.

    The writer side (annotation/tasks/calculate_sample_stats.calculate_cohort_stats)
    uses the same canonical_filter_key helper to populate filter-keyed rows. A
    silent mismatch between the two would produce cache misses forever.
"""
from typing import Iterable, Optional, Union

from django.core.exceptions import ObjectDoesNotExist

from analysis.models.enums import TrioInheritance
from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from library.utils.json_utils import canonical_filter_key
from snpdb.models import CohortGenotypeStats
from snpdb.models.models_enums import BuiltInFilters


class _Uncacheable:
    """ Sentinel returned by node handlers for configs we don't precompute. """


UNCACHEABLE = _Uncacheable()
FilterKey = Union[None, str, _Uncacheable]


LABEL_TO_MODEL = {
    BuiltInFilters.TOTAL: CohortGenotypeStats,
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: CohortGenotypeVariantAnnotationStats,
    BuiltInFilters.OMIM: CohortGenotypeGeneAnnotationStats,
    BuiltInFilters.CLINVAR: CohortGenotypeClinVarAnnotationStats,
}


class FilterKeyHandler:
    """ A node-class-specific helper that:
          1. Given a node, returns the filter_key (or UNCACHEABLE) describing its config.
          2. Given a cohort, returns the filter_keys the calc function should pre-populate.
    """

    def filter_key_for_node(self, node) -> FilterKey:
        raise NotImplementedError

    def filter_keys_to_precompute(self, cohort) -> Iterable[Optional[str]]:
        return ()


class NoFilterHandler(FilterKeyHandler):
    """ Node only ever wants the raw aggregate row (filter_key=NULL). """

    def filter_key_for_node(self, node) -> FilterKey:
        return None

    def filter_keys_to_precompute(self, cohort) -> Iterable[Optional[str]]:
        return [None]


class TrioInheritanceHandler(FilterKeyHandler):
    # 5 inheritance modes precomputed including compound het.
    CACHED_MODES = {
        TrioInheritance.DENOVO: "denovo",
        TrioInheritance.RECESSIVE: "autosomal_recessive",
        TrioInheritance.DOMINANT: "autosomal_dominant",
        TrioInheritance.XLINKED_RECESSIVE: "x_linked_recessive",
        TrioInheritance.COMPOUND_HET: "compound_het",
    }

    def filter_key_for_node(self, node) -> FilterKey:
        # TrioNode always has an inheritance set.
        mode = self.CACHED_MODES.get(TrioInheritance(node.inheritance))
        if mode is None:
            return UNCACHEABLE
        return canonical_filter_key({"inheritance": mode})

    def filter_keys_to_precompute(self, cohort) -> Iterable[Optional[str]]:
        if not _cohort_has_trio(cohort):
            return [None]
        return [None] + [canonical_filter_key({"inheritance": m}) for m in self.CACHED_MODES.values()]


def _cohort_has_trio(cohort) -> bool:
    return cohort.trio_set.exists()


def get_handler_for_node(node) -> FilterKeyHandler:
    # Local imports avoid eager imports of node modules during _stats_cache import.
    from analysis.models.nodes.sources.cohort_node import CohortNode
    from analysis.models.nodes.sources.pedigree_node import PedigreeNode
    from analysis.models.nodes.sources.sample_node import SampleNode
    from analysis.models.nodes.sources.trio_node import TrioNode

    handlers = {
        SampleNode: NoFilterHandler(),
        CohortNode: NoFilterHandler(),
        TrioNode: TrioInheritanceHandler(),
        PedigreeNode: NoFilterHandler(),
    }
    return handlers.get(type(node), NoFilterHandler())


def get_filter_keys_to_precompute_for_cohort(cohort) -> list[Optional[str]]:
    """ Used by the writer (calculate_cohort_stats) to know which filter_key
        buckets to populate beyond the raw aggregate (None) row. Trios get the
        5 inheritance keys; everything else gets just None. """
    keys: list[Optional[str]] = [None]
    if _cohort_has_trio(cohort):
        keys.extend(canonical_filter_key({"inheritance": m})
                    for m in TrioInheritanceHandler.CACHED_MODES.values())
    return keys


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
    from annotation.tasks.calculate_sample_stats import calculate_cohort_stats_task
    calculate_cohort_stats_task.apply_async(
        args=[cohort.pk, annotation_version.pk],
        queue="annotation_workers",
    )
