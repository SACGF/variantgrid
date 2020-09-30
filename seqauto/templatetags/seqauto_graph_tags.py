from django.template import Library
import json
import uuid

from seqauto.seqauto_stats import group_enrichment_kits_df


register = Library()


@register.inclusion_tag("seqauto/tags/sample_enrichment_kits_graph.html")
def sample_enrichment_kits_graph(sample_enrichment_kits_df, title, by_column):
    enrichment_kits_over_time, enrichment_kit_labels = group_enrichment_kits_df(sample_enrichment_kits_df, by_column, max_groups=10)

    return {'title': title,
            'uuid': uuid.uuid4(),
            'enrichment_kits_over_time': enrichment_kits_over_time,
            'enrichment_kit_labels': enrichment_kit_labels}
