import json

from django.template import Library

from analysis.models.nodes.analysis_node import NodeAlleleFrequencyFilter
from analysis.serializers import NodeAlleleFrequencyFilterSerializer

register = Library()


@register.inclusion_tag("analysis/tags/allele_frequency_tag.html")
def allele_frequency_controls(node):
    allele_frequency = get_allele_frequency_filter_json_for_node(node)

    return {'allele_frequency': allele_frequency}


def get_allele_frequency_filter_json_for_node(node):
    json_data = None
    try:
        serializer = NodeAlleleFrequencyFilterSerializer(node.nodeallelefrequencyfilter)
        json_data = json.dumps(serializer.data)
    except NodeAlleleFrequencyFilter.DoesNotExist:
        pass
    return json_data
