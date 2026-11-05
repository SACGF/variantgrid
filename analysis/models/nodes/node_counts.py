from typing import Optional

from django.db.models import Count
from django.db.models.query_utils import Q

from analysis.exceptions import NonFatalNodeError
from snpdb.enums import BuiltInFilters, TagFilter
from snpdb.models import Tag, Variant
from snpdb.models.models_user_settings import AbstractNodeCountSettings
from snpdb.utils import get_tag_styles_and_colors

# Tags the user hasn't given a colour still need a visible legend swatch
DEFAULT_TAG_NODE_COUNT_COLOR = "#888888"


def get_analysis_tagged_variant_ids(analysis, tag_ids: list[str]) -> list[int]:
    """ Variants carrying any of tag_ids in this analysis. Tagging is manual so this stays small -
        small enough to pass around as a list, @see TagNode.tagged_variants_q which does the same """
    # Go via Variant rather than VariantTag as analysis.models.models_variant_tag imports AnalysisNode,
    # which imports this module
    tagged_qs = Variant.objects.filter(varianttag__analysis=analysis, varianttag__tag__in=tag_ids)
    return list(tagged_qs.values_list("pk", flat=True).distinct())


def is_extra_filter(extra_filters) -> bool:
    """ Whether extra_filters narrows the node - "default" (and None) mean show the whole thing """
    return bool(TagFilter.get_tag_ids(extra_filters)) or extra_filters in dict(BuiltInFilters.FILTER_CHOICES)


def get_label_q(node, label) -> Q:
    """ A tag label's analysis-scoped tagged variants, otherwise the filter of the node the count label
        is configured as @see AnalysisNode.count_label_node """
    if tag_ids := TagFilter.get_tag_ids(label):
        return Q(pk__in=get_analysis_tagged_variant_ids(node.analysis, tag_ids))
    return node.count_label_q(node.analysis, label)


def get_node_extra_filters_q(node, extra_filters) -> Optional[Q]:
    """ The node's own filter for an 'extra_filters' selection, or None if it selects everything.
        A global TagNode's tags reach outside the analysis, so it decides its own tag scope -
        everything else is analysis-scoped, which is what the DAG's node counts are """
    if not is_extra_filter(extra_filters):
        return None
    if tag_ids := TagFilter.get_tag_ids(extra_filters):
        if (node_tag_q := node.get_extra_filters_tag_q(tag_ids)) is not None:
            return node_tag_q
    return get_label_q(node, extra_filters)


def get_extra_filters_count(node, extra_filters) -> Optional[int]:
    """ How many rows the node shows under extra_filters, or None if we can't say cheaply. There
        being no filter is one of those - the node's own count is what covers that """
    if not is_extra_filter(extra_filters):
        return None
    tag_ids = TagFilter.get_tag_ids(extra_filters)
    node_tag_q = node.get_extra_filters_tag_q(tag_ids) if tag_ids else None
    if node_tag_q is None:
        # Stored counts are analysis-scoped, so they only speak for an analysis-scoped filter
        if (count := node.node_version.counts.get(extra_filters)) is not None:
            return count
    if tag_ids:
        # Tags narrow the node to a handful of variants, so counting them exactly here is cheap
        q = node_tag_q if node_tag_q is not None else get_label_q(node, extra_filters)
        try:
            return node.get_queryset().filter(q).count()
        except NonFatalNodeError:
            pass  # An ancestor isn't ready - the count comes back when the node reloads
    return None


def get_node_count_colors(css_property):
    """ Returns a list of tuples with last element being a dict,
        css_property of "color" = [('ClinVar', {color: #ff0000}), etc] """

    node_count_colors = []
    for label, color in BuiltInFilters.COLORS:
        node_count_colors.append((label, {css_property: color}))

    return node_count_colors


def get_tag_node_count_colors(user, css_property, tag_config_collection=None):
    """ Tag node counts are drawn in the resolved tag config's colours - same shape as get_node_count_colors() """
    _, user_tag_colors = get_tag_styles_and_colors(user, tag_config_collection=tag_config_collection)
    return [(TagFilter.label(tag_id), {css_property: color or DEFAULT_TAG_NODE_COUNT_COLOR})
            for tag_id, color in user_tag_colors.items()]


def get_node_counts_mine_and_available(analysis):
    """ (mine, available built in filters, available tags) - mine keeps the user's order, the rest stay
        in their natural order. Tags are split out as there are far more of them than built in filters """
    node_count_types = analysis.get_node_count_types()

    my_choices = [x[0] for x in node_count_types]
    filter_choices = [x[0] for x in BuiltInFilters.CHOICES]
    tag_choices = [TagFilter.label(tag_id) for tag_id in Tag.live_qs().order_by("pk").values_list("pk", flat=True)]

    def _node_count_list(choices):
        node_counts_list = []
        for node_count in choices:
            if description := AbstractNodeCountSettings.get_node_count_description(node_count):
                node_counts_list.append({"pk": node_count,
                                         "css_classes": 'node-count-legend-' + node_count,
                                         "tag_id": TagFilter.get_tag_id(node_count),
                                         "description": description})
        return node_counts_list

    mine = set(my_choices)
    return (_node_count_list(my_choices),
            _node_count_list([c for c in filter_choices if c not in mine]),
            _node_count_list([c for c in tag_choices if c not in mine]))


def get_node_counts_and_labels_dict(node, counts_to_get):
    return _aggregate_node_counts(node, counts_to_get)


def get_tagged_variant_ids_by_label(analysis, tag_labels) -> dict[str, list[int]]:
    """ The tagged variants behind each tag node count label. The same for every node in the analysis,
        so look them up once and pass them to get_tag_node_counts_dict() per node """
    return {label: get_analysis_tagged_variant_ids(analysis, [TagFilter.get_tag_id(label)])
            for label in tag_labels}


def get_tag_node_counts_dict(node, tagged_variant_ids_by_label: dict[str, list[int]]) -> dict[str, int]:
    """ Tags change without bumping node versions, so these are recounted in place. Restricting the
        node queryset to the (small) tagged variants first makes this far cheaper than a node reload """
    all_tagged_variant_ids = set()
    for variant_ids in tagged_variant_ids_by_label.values():
        all_tagged_variant_ids.update(variant_ids)

    if not all_tagged_variant_ids:
        return dict.fromkeys(tagged_variant_ids_by_label, 0)

    qs = node.get_queryset().filter(pk__in=all_tagged_variant_ids)
    return _count_by_label(qs, {label: Q(pk__in=variant_ids)
                                for label, variant_ids in tagged_variant_ids_by_label.items()})


def _aggregate_node_counts(node, counts_to_get) -> dict[str, int]:
    q_by_label = {}
    for count_type in counts_to_get:
        if count_type == BuiltInFilters.TOTAL:
            q = None
        else:
            q = get_label_q(node, count_type)
        q_by_label[count_type] = q
    return _count_by_label(node.get_queryset(), q_by_label)


def _count_by_label(qs, q_by_label: dict[str, Optional[Q]]) -> dict[str, int]:
    """ Labels carry tag ids, which can hold characters Django refuses in a column alias (whitespace,
        quotes, ';') - so aggregate under generated aliases and map the counts back """
    labels_by_alias = {f"count_{i}": label for i, label in enumerate(q_by_label)}
    counts = qs.aggregate(**{alias: Count("pk", filter=q_by_label[label])
                             for alias, label in labels_by_alias.items()})
    return {label: counts[alias] for alias, label in labels_by_alias.items()}
