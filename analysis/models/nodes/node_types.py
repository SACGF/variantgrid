from collections import defaultdict
from django.db import models

from analysis.models.nodes.analysis_node import AnalysisNode
from library.utils import get_subclasses


class NodeGraphType(models.Model):
    name = models.TextField()
    required_column = models.TextField(default='id')
    graph_class = models.TextField()
    description = models.TextField()

    def __str__(self):
        return self.name


class NodeHelp(models.Model):
    node_class_name = models.TextField()
    help_text = models.TextField()

    def __str__(self):
        text = self.help_text
        if len(text) > 20:
            text = text[:20] + "..."
        return f"{self.node_class_name}/{text}"

    @staticmethod
    def get_node_help_dict():
        node_help_dict = {}

        # Set to empty in case it's not in DB
        for node_class_name in get_node_types_hash():
            node_help_dict[node_class_name] = ''

        # Get everything out of DB at once
        node_help = NodeHelp.objects.all().values_list("node_class_name", "help_text")
        for node_class_name, help_text in node_help:
            node_help_dict[node_class_name] = help_text

        return node_help_dict


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


def get_nodes_by_classification():
    node_types = get_node_types_hash()
    nodes = defaultdict(list)

    node_help_dict = NodeHelp.get_node_help_dict()

    for node_class_name, node_class in node_types.items():
        node = node_class()
        node_class = node.get_class_name()

        max_inputs = node.max_inputs
        if max_inputs == -1:
            max_inputs = 'Unlimited'

        data = {"node_class": node_class,
                "node_class_name": node_class_name,
                "help_text": node_help_dict[node_class_name],
                "min_inputs": node.min_inputs,
                "max_inputs": max_inputs}
        classification = node.get_node_classification()
        nodes[classification].append(data)

    return nodes
