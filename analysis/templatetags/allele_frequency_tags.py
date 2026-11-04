from django.conf import settings
from django.template import Library

from analysis.models.nodes.node_allele_frequency import NodeAlleleFrequencyFilter
from analysis.serializers import NodeAlleleFrequencyFilterSerializer

register = Library()


@register.inclusion_tag("analysis/tags/allele_frequency_tag.html")
def allele_frequency_controls(node):
    allele_frequency = get_allele_frequency_filter_data_for_node(node)

    return {
        "VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT": settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT,
        'allele_frequency': allele_frequency
    }


def get_allele_frequency_filter_data_for_node(node):
    try:
        return NodeAlleleFrequencyFilterSerializer(node.nodeallelefrequencyfilter).data
    except NodeAlleleFrequencyFilter.DoesNotExist:
        return None
