from collections import defaultdict
from typing import Dict, List

from django.db import models

from analysis.models.nodes.analysis_node import AnalysisNode
from library.utils import get_subclasses


class NodeGraphType(models.Model):
    name = models.TextField()
    required_column = models.TextField(default='id')
    graph_class = models.TextField()
    description = models.TextField()

    def __str__(self) -> str:
        return self.name


def get_node_types_hash():
    """ hash by get_node_class_label() """
    subclasses = {}
    for sc in get_subclasses(AnalysisNode):
        if sc.disabled or sc._meta.abstract:
            continue

        label = sc.get_node_class_label()
        subclasses[label] = sc
    return subclasses


def get_node_types_hash_by_class_name():
    node_types = get_node_types_hash()

    node_types_hash_by_class = {}
    for node_class in node_types.values():
        class_name = node_class().__class__.__name__
        node_types_hash_by_class[class_name] = node_class
    return node_types_hash_by_class


def get_nodes_by_classification() -> Dict[str, List]:
    node_types = get_node_types_hash()
    nodes = defaultdict(list)

    for node_class_name, node_class in node_types.items():
        node = node_class()
        classification = node.get_node_classification()
        nodes[classification].append(node.get_class_name())

    return nodes
