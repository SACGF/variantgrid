from collections import defaultdict
from dataclasses import asdict

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
    return {node_class.__name__: node_class for node_class in get_node_types_hash().values()}


def get_nodes_by_classification() -> dict[str, list]:
    """ Add node dropdown rows - one per node class """
    nodes = defaultdict(list)
    for class_label, node_class in get_node_types_hash().items():
        node = node_class()
        classification = node.get_node_classification()
        nodes[classification].append({
            "class_name": node.get_class_name(),  # What node_create is passed, and the accent colour
            "class_label": class_label,
            "classification": classification,  # add node dropdown colours icons like the cards
            "icon": asdict(node_class.get_node_class_icon()),
        })
    return nodes


def get_node_display_data_by_class_name() -> dict[str, dict]:
    """ Icons/labels for the add node dropdown, keyed by the <select> option values """
    return {data["class_name"]: data
            for nodes in get_nodes_by_classification().values()
            for data in nodes}
