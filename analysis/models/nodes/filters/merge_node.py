import logging
import operator
from functools import reduce
from typing import Optional, Iterable

from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from snpdb.models import VariantCollectionRecord, lazy


class MergeNode(AnalysisNode):
    min_inputs = 1
    max_inputs = AnalysisNode.PARENT_CAP_NOT_SET

    def __init__(self, *args, **kwargs):
        kwargs["parents_should_cache"] = True
        super().__init__(*args, **kwargs)

    def modifies_parents(self):
        return self._num_unique_parents_in_queryset > 1

    def _queryset_requires_distinct(self):
        # Need distinct when we have multiple Q's OR'd
        # No need if there are multiple parents but all cached
        return len(self._unique_parent_q_set) > 1

    def get_cache_task_args_objs_set(self, force_cache=False):
        cache_task_args_objs = set()
        if self.is_valid():
            for p in self.get_non_empty_parents(require_parents_ready=False):
                # Ensure that parents all have cache. This is only true if the cache is READY - it may be being
                # written in another task.
                if p.node_cache is None:
                    logging.warning("Parent %s didn't have cache!!!", p)
                    cache_task_args_objs.update(p.get_cache_task_args_objs_set(force_cache=True))
        return cache_task_args_objs

    @lazy
    def _num_unique_parents_in_queryset(self):
        return len(set(p.get_q() for p in self.get_non_empty_parents(require_parents_ready=False)))

    def get_single_parent(self):
        """ Override so we can use get_grid_node_id_and_version
            The query AND the input samples must be the same """

        if self._num_unique_parents_in_queryset == 1:
            my_sample_ids = self.get_sample_ids()
            for parent in self.get_non_empty_parents():  # As only 1 unique, can take 1st we find
                if parent.get_sample_ids() == my_sample_ids:
                    # print(f"MergeNode using {parent} as single unmodified parent!")
                    return parent
        return super().get_single_parent()  # Will throw exception due to multiple samples

    @lazy
    def _unique_parent_q_set(self) -> Iterable[Q]:
        """ Set to not duplicate Q's (eg from node that doesn't modify parents and that node's parent) """
        parent_variant_collection_ids = set()
        unique_parent_q_set = set()
        for parent in self.get_non_empty_parents():
            if parent.node_cache:
                parent_variant_collection_ids.add(parent.node_cache.variant_collection_id)
            else:
                unique_parent_q_set.add(parent.get_q())

        if parent_variant_collection_ids:
            variants_from_cache_qs = VariantCollectionRecord.objects.filter(
                variant_collection__in=sorted(parent_variant_collection_ids))  # Sort so Q object strings hash the same
            unique_parent_q_set.add(Q(pk__in=variants_from_cache_qs.values_list("variant_id")))

        # Q objects for the same ultimate query don't hash the same, so remove duplicates via string equality
        # This can save a distinct
        if len(unique_parent_q_set) > 1:
            q_by_string = {str(q): q for q in unique_parent_q_set}
            return q_by_string.values()
        return unique_parent_q_set

    def get_parent_arg_q_dict(self):
        if self._unique_parent_q_set:
            q = reduce(operator.or_, self._unique_parent_q_set)
        else:
            q = self.q_none()
        return q

    def _get_node_q(self) -> Optional[Q]:
        return None  # No extra filtering needed after get_parent_q()

    def _get_method_summary(self):
        parent_names = ','.join([p.name for p in self.get_parent_subclasses()])
        return f"Merged from parents: {parent_names}"

    def get_node_name(self):
        return "Merge"

    @staticmethod
    def get_help_text() -> str:
        return "Merge variants from multiple parents"

    @staticmethod
    def get_node_class_label():
        return "Merge"
