"""
    Reader-side cache lookups for the CohortGenotype*Stats family + node-class
    handlers that translate node configs to canonical filter_keys.

    The writer side (annotation/tasks/calculate_sample_stats.calculate_cohort_stats)
    uses the same canonical_filter_key helper to populate filter-keyed rows. A
    silent mismatch between the two would produce cache misses forever.
"""
from collections.abc import Iterable
from typing import Optional, Union

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


class SampleNodeHandler(NoFilterHandler):
    """ A group level SampleNode spans genotype collections, so there's no single cohort's stats
        row to read - its counts are live queries. """

    def filter_key_for_node(self, node) -> FilterKey:
        if node.is_group_level:
            return UNCACHEABLE
        return None


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
        # The precomputed buckets require both parents to have a zygosity call (@see _trio_predicates).
        # require_zygosity=False also matches parent no-calls, which no bucket represents.
        if not node.require_zygosity:
            return UNCACHEABLE
        # TrioNode always has an inheritance set.
        mode = self.CACHED_MODES.get(TrioInheritance(node.inheritance))
        if mode is None:
            return UNCACHEABLE
        return inheritance_filter_key(mode, node.trio)

    def filter_keys_to_precompute(self, cohort) -> Iterable[Optional[str]]:
        trio = _cohort_trio(cohort)
        if trio is None:
            return [None]
        return [None] + [inheritance_filter_key(m, trio) for m in self.CACHED_MODES.values()]


class DuoInheritanceHandler(FilterKeyHandler):
    """ Nothing is precomputed for duos - the buckets are built off a cohort's Trio (@see
        _trio_predicates), so every DuoNode counts live rather than reading the raw cohort row,
        which knows nothing about its inheritance mode. """

    def filter_key_for_node(self, node) -> FilterKey:
        return UNCACHEABLE


class QuadInheritanceHandler(FilterKeyHandler):
    """ Nothing is precomputed for quads - the buckets are built off a cohort's Trio (@see
        _trio_predicates), so every QuadNode counts live rather than reading the raw cohort row,
        which knows nothing about its inheritance mode. """

    def filter_key_for_node(self, node) -> FilterKey:
        return UNCACHEABLE


def inheritance_filter_key(mode: str, trio) -> str:
    """ The key a mode's bucket is stored under - the writer and the readers both come through here
        so they can't drift apart. autosomal_dominant's predicate branches on which parents are
        affected (@see _trio_predicates), so those flags belong in its key: editing them misses the
        old bucket and recomputes, rather than reading one built for the previous flags. """
    key = {"inheritance": mode}
    if mode == "autosomal_dominant":
        key["mother_affected"] = bool(trio.mother_affected)
        key["father_affected"] = bool(trio.father_affected)
    return canonical_filter_key(key)


def _cohort_trio(cohort):
    return cohort.trio_set.first()


def get_handler_for_node(node) -> FilterKeyHandler:
    # Local imports avoid eager imports of node modules during stats_cache import.
    from analysis.models.nodes.sources.cohort_node import CohortNode
    from analysis.models.nodes.sources.duo_node import DuoNode
    from analysis.models.nodes.sources.pedigree_node import PedigreeNode
    from analysis.models.nodes.sources.quad_node import QuadNode
    from analysis.models.nodes.sources.sample_node import SampleNode
    from analysis.models.nodes.sources.trio_node import TrioNode

    handlers = {
        SampleNode: SampleNodeHandler(),
        CohortNode: NoFilterHandler(),
        TrioNode: TrioInheritanceHandler(),
        DuoNode: DuoInheritanceHandler(),
        QuadNode: QuadInheritanceHandler(),
        PedigreeNode: NoFilterHandler(),
    }
    return handlers.get(type(node), NoFilterHandler())


def get_filter_keys_to_precompute_for_cohort(cohort) -> list[Optional[str]]:
    """ Used by the writer (calculate_cohort_stats) to know which filter_key
        buckets to populate beyond the raw aggregate (None) row. Trios get the
        5 inheritance keys; everything else gets just None. """
    return list(TrioInheritanceHandler().filter_keys_to_precompute(cohort))


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
    from annotation.tasks.calculate_sample_stats import enqueue_cohort_stats_recompute
    enqueue_cohort_stats_recompute(cohort, annotation_version)
