import operator
from collections import defaultdict, Counter
from functools import cached_property, reduce
from typing import Optional

from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from variantgrid import settings


class MergeNode(AnalysisNode):
    min_inputs = 1
    max_inputs = AnalysisNode.PARENT_CAP_NOT_SET

    def modifies_parents(self):
        return self._num_unique_parents_in_queryset > 1

    @cached_property
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

    @staticmethod
    def _split_common_filters(parent_arg_q_dict):
        """ Find common filters - this only works on arg = None (as those can be re-composed) """

        arg_q_nodes = defaultdict(set)
        all_q_by_hash = {}

        for parent, arg_q_dict in parent_arg_q_dict.items():
            # We only use arg=None as that can be combined
            for q_hash, q in arg_q_dict.get(None, {}).items():
                arg_q_nodes[q_hash].add(parent)
                all_q_by_hash[q_hash] = q

        or_list = []
        non_combine_parents = set(parent_arg_q_dict.keys())

        if arg_q_nodes:
            # Find ones that are common
            most_common = sorted(arg_q_nodes.items(), key=lambda x: len(x[1]), reverse=True)[0]
            combine_q_hash, combine_parents = most_common
            if len(combine_parents) > 1:
                non_combine_parents -= combine_parents
                combine_parent_arg_q_dict = {p: parent_arg_q_dict[p] for p in combine_parents}

                extract_arg_q_hash = {None: {combine_q_hash}}
                MergeNode._remove_arg_q_hash(combine_parent_arg_q_dict, extract_arg_q_hash)
                combine_q = all_q_by_hash[combine_q_hash]
                if common_q := MergeNode._split_common_filters(combine_parent_arg_q_dict):
                    combine_q &= common_q
                or_list.append(combine_q)

        for parent in non_combine_parents:
            arg_q_dict = parent_arg_q_dict[parent]
            # If there is something else other than None - then we need to run the full queryset
            non_none_keys = [k for k in arg_q_dict.keys() if k is not None]

            if non_none_keys:
                print(f"{non_none_keys=}")
                # We don't pass in arg_q_dict (ie run all where clauses in inner query)
                # This has worse best-case performance but better worse case performance
                qs = parent.get_queryset(disable_cache=True)
                variant_ids = qs.values_list("pk", flat=True)
                or_list.append(Q(pk__in=variant_ids))
            else:
                if remaining_q_dict := arg_q_dict.get(None):
                    merged_q = reduce(operator.and_, remaining_q_dict.values())
                    or_list.append(merged_q)

        q = None
        if or_list:
            q = reduce(operator.or_, or_list)
        return q

    @staticmethod
    def _remove_arg_q_hash(parent_arg_q_dict, extract_arg_q_hash):
        for arg_q_dict in parent_arg_q_dict.values():
            for arg, q_hash_set in extract_arg_q_hash.items():
                if q_dict := arg_q_dict.get(arg, {}):
                    for q_hash in q_hash_set:
                        q_dict.pop(q_hash, None)
                    if not q_dict:
                        del arg_q_dict[arg]

    def _get_merged_q_dict(self, parent_arg_q_dict):
        all_q_by_hash = {}
        filtered_relation_count = defaultdict(Counter)
        for arg_q_dict in parent_arg_q_dict.values():
            for arg, q_dict in arg_q_dict.items():
                for q_hash, q in q_dict.items():
                    all_q_by_hash[q_hash] = q
                    if arg is not None:
                        filtered_relation_count[arg][q_hash] += 1

        num_parents = len(parent_arg_q_dict)
        extract_arg_q_hash = defaultdict(set)
        for arg, q_hash_count in filtered_relation_count.items():
            for q_hash, count in q_hash_count.items():
                if count == num_parents:
                    extract_arg_q_hash[arg].add(q_hash)

        arg_q_dict = {}
        for arg, q_hash_set in extract_arg_q_hash.items():
            arg_q_dict[arg] = {q_hash: all_q_by_hash[q_hash] for q_hash in q_hash_set}

        self._remove_arg_q_hash(parent_arg_q_dict, extract_arg_q_hash)

        if q := self._split_common_filters(parent_arg_q_dict):
            arg_q_dict[None] = {self._get_node_q_hash(): q}
        return arg_q_dict

    def _get_arg_q_dict_from_parents_and_node(self):
        parent_arg_q_dict = {}
        for parent in self.get_non_empty_parents():
            if settings.ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX and \
                    parent.count <= settings.ANALYSIS_NODE_MERGE_STORE_ID_SIZE_MAX:
                variant_ids = list(parent.get_queryset().values_list("pk", flat=True))
                q = Q(pk__in=variant_ids)
                arg_q_dict = {None: {q: q}}
            else:
                arg_q_dict = parent.get_arg_q_dict(disable_cache=True)
            parent_arg_q_dict[parent] = arg_q_dict
        return self._get_merged_q_dict(parent_arg_q_dict)

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
