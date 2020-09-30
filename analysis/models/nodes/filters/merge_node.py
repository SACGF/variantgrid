from functools import reduce
import operator
from typing import Optional

from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode


class MergeNode(AnalysisNode):
    min_inputs = 1
    max_inputs = AnalysisNode.PARENT_CAP_NOT_SET

    def __init__(self, *args, **kwargs):
        kwargs["parents_should_cache"] = True
        super().__init__(*args, **kwargs)

    def modifies_parents(self):
        return self.num_parents > 1

    def _queryset_requires_distinct(self):
        return self.modifies_parents()

    def get_cache_task_args_objs_set(self, force_cache=False):
        cache_task_args_objs = set()
        if self.is_valid():
            for p in self.get_non_empty_parents(require_parents_ready=False):
                # Ensure that parents all have cache. This is only true if the cache is READY - it may be being
                # written in another task.
                if p.node_cache is None:
                    print(f"Parent {p} didn't have cache!!!")
                    cache_task_args_objs.update(p.get_cache_task_args_objs_set(force_cache=True))
        return cache_task_args_objs

    @property
    def num_parents(self):
        """ Num effective parents (who contribute to query) """
        return len(self.get_unique_parent_q())

    def get_single_parent(self):
        """ Override so we can use get_grid_node_id_and_version
            The query AND the input samples must be the same """

        if self.num_parents == 1:
            my_sample_ids = self.get_sample_ids()
            for parent in self.get_non_empty_parents():
                if parent.get_sample_ids() == my_sample_ids:
                    # print(f"MergeNode using {parent} as single unmodified parent!")
                    return parent

        return super().get_single_parent()  # Will throw exception due to multiple samples

    def get_unique_parent_q(self):
        return set(p.get_q() for p in self.get_non_empty_parents())

    def get_parent_q(self):
        parent_q = self.get_unique_parent_q()
        if parent_q:
            # Will call distinct() on this if more than 1 parent_q (see above)
            q = reduce(operator.or_, parent_q)
        else:
            q = self.q_none()
        return q

    def _get_node_q(self) -> Optional[Q]:
        return None  # All handled in get_parent_q()

    def _get_method_summary(self):
        parent_names = ','.join([p.name for p in self.get_parent_subclasses()])
        return f"Merged from parents: {parent_names}"

    def get_node_name(self):
        return "Merge"

    @staticmethod
    def get_node_class_label():
        return "Merge"
