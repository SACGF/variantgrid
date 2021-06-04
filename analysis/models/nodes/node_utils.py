from celery.canvas import Signature
from collections import defaultdict
from django.db.models.query_utils import Q
from toposort import toposort
import random

from analysis.models.nodes.analysis_node import AnalysisNode, AnalysisEdge


def get_nodes_by_id(nodes_qs):
    nodes_by_id = {}

    for node in nodes_qs:
        nodes_by_id[node.id] = node
    return nodes_by_id


def get_parent_value_dag_dictionary(nodes):
    """ Create dict that's used for toposort:

        Dependencies are expressed as a dictionary whose keys are items
        and whose values are a set of dependent items aka 'parent value DAG' """

    graph = defaultdict(set)

    for node_id in nodes:
        graph[node_id] = set()

    query = Q(parent__in=nodes) | Q(child__in=nodes)
    qs = AnalysisEdge.objects.filter(query).values_list('parent_id', 'child_id')  # @UndefinedVariable
    for parent_id, child_id in qs:
        graph[child_id].add(parent_id)
    return graph


def get_toposorted_nodes(nodes_qs):
    nodes = get_nodes_by_id(nodes_qs)
    parent_value_data = get_parent_value_dag_dictionary(nodes)
    return get_toposorted_nodes_from_parent_value_data(nodes, parent_value_data)


def get_toposorted_nodes_from_parent_value_data(nodes, parent_value_data):
    topo_sorted_nodes = []
    for grp in toposort(parent_value_data):
        nodes_group = []
        for node_id in grp:
            node = nodes.get(node_id)
            if node:
                nodes_group.append(node)
        if nodes_group:
            topo_sorted_nodes.append(nodes_group)

    return topo_sorted_nodes


def update_analysis(analysis_id):
    """ Launches async job to update analysis """
    task = Signature("analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks", args=(analysis_id,))
    task.apply_async()


def reload_analysis_nodes(analysis_id):
    for node in AnalysisNode.objects.filter(analysis_id=analysis_id).select_subclasses():
        node.update_children = False  # Will get them all in loop
        node.appearance_dirty = True
        node.queryset_dirty = True
        node.save()
    return update_analysis(analysis_id)


def get_rendering_dict(node):
    node_class = node.get_class_name()
    css_classes = ["node-overlay"] + node.get_css_classes()  # Needs node-overlay to be able to find it
    node_args = node.get_rendering_args()

    if node.pk:
        node_id = int(node.pk)
    else:
        node_id = None

    style = f"left: {node.x}px; top: {node.y}px"
    attributes = {
        "node_id": node_id,
        "node_class": node.get_node_class_label(),
        "version_id": node.version,
        "appearance_version_id": node.appearance_version,
        "id": node.get_css_id(),
        "style": style,
        "input_endpoint": node.max_inputs != 0,
        "output_endpoint": True,  # Can always have output
        "x": node.x,
        "y": node.y
    }
    return {
        "attributes": attributes,
        "node_class": node_class,
        "overlay_css_classes": " ".join(css_classes),
        "name": node.name,
        "args": node_args
    }


def create_node(node_class, analysis, **kwargs):
    save = kwargs.pop("save", True)
    node_kwargs = {"x": 10 + random.random() * 50,
                   "y": 10 + random.random() * 20}
    node_kwargs.update(kwargs)
    node_kwargs["analysis"] = analysis
    if save:
        node = node_class.objects.create(**node_kwargs)
    else:
        node = node_class(**node_kwargs)
    return node
