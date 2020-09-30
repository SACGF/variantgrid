from collections import defaultdict
from django import template
import json
import logging

from analysis.models import NodeVariant
from analysis.models.nodes.node_utils import get_rendering_dict
from library import tag_utils

register = template.Library()


@register.filter(name='get_class')
def get_class(value):
    return value.get_class_name()


class AnalysisNodesTemplateNode(template.Node):

    def __init__(self, nodes):
        self.nodes_variable = template.Variable(nodes)

    def render(self, context):
        try:
            nodes = self.nodes_variable.resolve(context)
            nodes_array = []
            for node in nodes:
                try:
                    node_rendering_dict = get_rendering_dict(node)
                    nodes_array.append(node_rendering_dict)
                except Exception as e:
                    logging.error(e)

            return json.dumps(nodes_array)
        except template.VariableDoesNotExist:
            return ''


class AnalysisNodeConnectionsJSNode(template.Node):
    """ Renders javascript to connect nodes in JSPlumb """

    def __init__(self, nodes):
        self.nodes_variable = template.Variable(nodes)

    def render(self, context):
        try:
            nodes = self.nodes_variable.resolve(context)
            connections = []
            for node in nodes:
                # NOTE: I am not sure why I have to call analysisnode_ptr but I do...
                if hasattr(node, "analysisnode_ptr"):
                    for parent in node.analysisnode_ptr.parents():
                        values = node.get_connection_data(parent)
                        connections.append(values)
            return json.dumps(connections)
        except template.VariableDoesNotExist:
            return ''


class NodeVariantNode(template.Node):

    def __init__(self, nodes):
        self.variable = template.Variable(nodes)

    def render(self, context):
        analysis = self.variable.resolve(context)

        node_variant = defaultdict(dict)
        node_variant_qs = NodeVariant.objects.filter(node__analysis=analysis).values_list('node_id', 'variant_id')
        for (node_id, variant_id) in node_variant_qs:
            node_variant[node_id][variant_id] = 1
        return json.dumps(node_variant)


@register.tag
def render_nodes_array(_parser, token):
    return AnalysisNodesTemplateNode(tag_utils.get_passed_object(token))


@register.tag
def render_nodes_connections_array(_parser, token):
    return AnalysisNodeConnectionsJSNode(tag_utils.get_passed_object(token))


@register.tag
def render_node_variant_dict(_parser, token):
    return NodeVariantNode(tag_utils.get_passed_object(token))
