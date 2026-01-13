import uuid

from django.template import Library

from seqauto.seqauto_stats import group_enrichment_kits_df

register = Library()


@register.inclusion_tag("seqauto/tags/sample_enrichment_kits_graph.html")
def sample_enrichment_kits_graph(sample_enrichment_kits_df, title, by_column, max_years: int = None):
    enrichment_kits_over_time, enrichment_kit_labels = group_enrichment_kits_df(sample_enrichment_kits_df, by_column,
                                                                                max_groups=10, max_years=max_years)

    if max_years is not None:
        title += f" (last {max_years} years)"

    return {
        'title': title,
        'uuid': uuid.uuid4(),
        'enrichment_kits_over_time': enrichment_kits_over_time,
        'enrichment_kit_labels': enrichment_kit_labels
    }
