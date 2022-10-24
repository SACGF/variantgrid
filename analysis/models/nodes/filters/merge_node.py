import operator
import time
from collections import defaultdict, Counter
from functools import reduce
from typing import Optional

from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from snpdb.models import lazy


class MergeNode(AnalysisNode):
    min_inputs = 1
    max_inputs = AnalysisNode.PARENT_CAP_NOT_SET

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache_node_q = False  # can't pickle querysets

    def modifies_parents(self):
        return self._num_unique_parents_in_queryset > 1

    @lazy
    def _num_unique_parents_in_queryset(self):
        parent_q_dicts = set()
        for p in self.get_non_empty_parents(require_parents_ready=False):
            key = tuple(((k, tuple(v)) for k, v in p.get_arg_q_dict().items()))
            parent_q_dicts.add(key)
        return len(parent_q_dicts)

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

    def _get_merged_q_dict(self, parent_arg_q_dict):
        """ This is the old way - used to compare against """
        q_or = []
        for parent, arg_q_dict in parent_arg_q_dict.items():
            qs = parent.get_queryset(disable_cache=True)
            q_or.append(Q(pk__in=qs.values_list("pk", flat=True)))

        return {
            None: {self._get_node_q_hash(): reduce(operator.or_, q_or)}
        }

    @staticmethod
    def _split_common_filters(parent_arg_q_dict):
        """ Find common filters """
        # Can you do it one node at a time?

        # find the node most frequently used by the parents


        # Find ones that are common
        # value = set of parent_id that shares Q
        arg_q_nodes = defaultdict(lambda: defaultdict(set))
        all_q_by_hash = {}

        for parent, arg_q_dict in parent_arg_q_dict.items():
            for arg, q_dict in arg_q_dict.items():
                for q_hash, q in q_dict.items():
                    arg_q_nodes[arg][q_hash].add(parent.pk)
                    all_q_by_hash[q_hash] = q

        # Find ones that are common
        parents_common = Counter()
        for values in arg_q_nodes.values():
            for parent_set in values.values():
                if len(parent_set) > 1:
                    parent_hash = tuple(sorted(parent_set))
                    parents_common[parent_hash] += 1

        if parents_common:
            parents, _count = sorted(parents_common.items(), key=operator.itemgetter(1), reverse=True)[0]
            parents = set(parents)
            common_q_list = []
            for parent, arg_q_dict in parent_arg_q_dict.items():
                if parent in parents:
                    pass # Put into filter, also remove it from parent list



            # Go through and extract out what we can combine vs what we can't
            # Put into an OR
            common = []
            q_or = []
            q = Q()
        else:
            # Nothing in common
            q_or = []
            for parent, arg_q_dict in parent_arg_q_dict.items():
                # We should also try just not passing in arg_q_dict here - in effect duplicating the filters
                # That are already applied on the outside - not sure what's better
                qs = parent.get_queryset(arg_q_dict=arg_q_dict, disable_cache=True)
                q_or.append(Q(pk__in=qs.values_list("pk", flat=True)))
            q = reduce(operator.or_, q_or)

        return q


    def _get_merged_q_dict2(self, parent_arg_q_dict):
        # value = set of parent_id that shares Q
        arg_q_nodes = defaultdict(lambda: defaultdict(set))
        all_q_by_hash = {}

        for parent, arg_q_dict in parent_arg_q_dict.items():
            for arg, q_dict in arg_q_dict.items():
                for q_hash, q in q_dict.items():
                    arg_q_nodes[arg][q_hash].add(parent.pk)
                    all_q_by_hash[q_hash] = q

        # Find ones that are common
        parents_common = Counter()
        for values in arg_q_nodes.values():
            for parent_set in values.values():
                if len(parent_set) > 1:
                    parent_hash = tuple(sorted(parent_set))
                    parents_common[parent_hash] += 1

        # Get the top one
        # divide into those that have it and don't
        # Those that don't - do the leaf query

        # Find the q_hash that has the highest count
        # divide the parents into those that have it and don't
        # Make an OR query with (   |  )
        # Add the query to
        # build an OR
        # extract it out
        #

        arg_q_dict = {}

        # arg_q_nodes
        # Maybe go and sort via count - take the highest
        for arg, q_dict in arg_q_nodes.items():
            print(arg)
            for q_hash, parents in q_dict.items():
                print(f"{q_hash}: {parents}")

        # TODO: Go and knock them out of the others, and leave only what's unique
        q_or = []
        for parent, arg_q_dict in parent_arg_q_dict.items():
            qs = parent.get_queryset(disable_cache=True)  # TODO: Pass in modified/unique arg_q_dict here
            q_or.append(Q(pk__in=qs.values_list("pk", flat=True)))

        arg_q_dict[None] = {self._get_node_q_hash(): reduce(operator.or_, q_or)}
        return arg_q_dict

    def _get_arg_q_dict_from_parents_and_node(self):
        # Go through and get the common things to all parents
        start = time.time()
        parent_arg_q_dict = {}
        for parent in self.get_non_empty_parents():
            arg_q_dict = parent.get_arg_q_dict(disable_cache=True)
            parent_arg_q_dict[parent] = arg_q_dict

        arg_q_dict = self._get_merged_q_dict(parent_arg_q_dict)
        # arg_q_dict = self._get_merged_q_dict(parent_arg_q_dict)

        end = time.time()
        print(f"merge calculations took {end-start} secs")
        return arg_q_dict

    def _get_node_q(self) -> Optional[Q]:
        raise NotImplementedError("This should never be called")

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
