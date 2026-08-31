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


def get_nodes_by_classification() -> dict[str, list]:
    """ Add node dropdown rows. A class contributes one row unless it declares several menu entries
        (SampleNode, one per source level) """
    nodes = defaultdict(list)

    for node_class in get_node_types_hash().values():
        node = node_class()
        classification = node.get_node_classification()
        for entry in node_class.get_menu_entries():
            nodes[classification].append({
                "key": entry.key,  # What node_create is passed
                "class_name": node.get_class_name(),  # Picks up the node's accent colour
                "class_label": entry.label,
                "class_label_short": node_class.get_node_class_label_short(),
                "classification": classification,  # add node dropdown colours icons like the cards
                "icon": asdict(entry.icon),
            })

    return nodes


def get_menu_entries_by_key() -> dict[str, tuple]:
    """ (node class, initial kwargs) for each add node dropdown row """
    entries = {}
    for node_class in get_node_types_hash().values():
        for entry in node_class.get_menu_entries():
            entries[entry.key] = (node_class, entry.initial_kwargs)
    return entries


def get_node_display_data_by_menu_key() -> dict[str, dict]:
    """ Icons/labels for the add node dropdown, keyed by the <select> option values """
    return {data["key"]: data
            for nodes in get_nodes_by_classification().values()
            for data in nodes}
