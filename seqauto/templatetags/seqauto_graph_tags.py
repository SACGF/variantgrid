import uuid
from typing import Optional

from django.template import Library

from seqauto.seqauto_stats import (
    VARIANTS_TYPE_SERIES_COL,
    VARIANTS_TYPE_SERIES_COLORS,
    group_enrichment_kits_df,
)

register = Library()


def _stacked_bar_context(groups, title, series_colors=None):
    return {
        'title': title,
        'uuid': uuid.uuid4(),
        'enrichment_kits_over_time': groups.data,
        'enrichment_kit_labels': groups.labels,
        'series_colors': series_colors or {},
        'collapsed_help': groups.collapsed_help,
    }


def _title_with_years(title, max_years):
    if max_years is not None:
        title += f" (last {max_years} years)"
    return title


@register.inclusion_tag("seqauto/tags/sample_enrichment_kits_graph.html")
def sample_enrichment_kits_graph(sample_enrichment_kits_df, title, by_column, max_years: Optional[int] = None):
    groups = group_enrichment_kits_df(sample_enrichment_kits_df, by_column, max_groups=10, max_years=max_years)
    return _stacked_bar_context(groups, _title_with_years(title, max_years))


@register.inclusion_tag("seqauto/tags/sample_enrichment_kits_graph.html")
def sample_variants_type_graph(sample_enrichment_kits_df, title, by_column, max_years: Optional[int] = None):
    """ Germline / somatic / unknown series from each sample's EnrichmentKit.sample_variants_type """
    groups = group_enrichment_kits_df(sample_enrichment_kits_df, by_column, max_years=max_years,
                                      group_column=VARIANTS_TYPE_SERIES_COL)
    return _stacked_bar_context(groups, _title_with_years(title, max_years), VARIANTS_TYPE_SERIES_COLORS)
