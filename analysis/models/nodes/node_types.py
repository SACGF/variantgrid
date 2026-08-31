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


def _walk_menu_entries():
    """ (node class, node instance, menu entry) for every add node dropdown row """
    for node_class in get_node_types_hash().values():
        node = node_class()
        for entry in node_class.get_menu_entries():
            yield node_class, node, entry


def get_nodes_by_classification() -> dict[str, list]:
    """ Add node dropdown rows. A class contributes one row unless it declares several menu entries
        (SampleNode, one per source level) """
    nodes = defaultdict(list)
    for _node_class, node, entry in _walk_menu_entries():
        classification = node.get_node_classification()
        nodes[classification].append({
            "key": entry.key,  # What node_create is passed
            "class_name": node.get_class_name(),  # Picks up the node's accent colour
            "class_label": entry.label,
            "classification": classification,  # add node dropdown colours icons like the cards
            "icon": asdict(entry.icon),
        })
    return nodes


def get_menu_entries_by_key() -> dict[str, tuple]:
    """ (node class, menu entry) for each add node dropdown row """
    return {entry.key: (node_class, entry) for node_class, _node, entry in _walk_menu_entries()}


def get_node_display_data_by_menu_key() -> dict[str, dict]:
    """ Icons/labels for the add node dropdown, keyed by the <select> option values """
    return {data["key"]: data
            for nodes in get_nodes_by_classification().values()
            for data in nodes}
