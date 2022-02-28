import logging

from django.contrib.auth.models import User
from django.db.models import Count
from django.db.models.query_utils import Q

from annotation.models.damage_enums import PathogenicityImpact
from classification.enums import ClinicalSignificance
from classification.models import Classification, GenomeBuild
from snpdb.models.models_enums import BuiltInFilters


def get_extra_filters_q(user: User, genome_build: GenomeBuild, extra_filters):
    if extra_filters == BuiltInFilters.CLINVAR:
        q = Q(clinvar__highest_pathogenicity__gte=4)
    elif extra_filters == BuiltInFilters.OMIM:
        q = Q(variantannotation__gene__geneannotation__omim_terms__isnull=False)
    elif extra_filters in [BuiltInFilters.CLASSIFIED, BuiltInFilters.CLASSIFIED_PATHOGENIC]:
        clinical_significance_list = None
        if extra_filters == BuiltInFilters.CLASSIFIED_PATHOGENIC:
            clinical_significance_list = [ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC]
        q = Classification.get_variant_q(user, genome_build, clinical_significance_list)
    elif extra_filters == BuiltInFilters.IMPACT_HIGH_OR_MODERATE:
        q = Q(variantannotation__impact__in=(PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE))
    elif extra_filters == BuiltInFilters.COSMIC:
        q = Q(variantannotation__cosmic_id__isnull=False)
    else:
        logging.warning("get_extra_filters_q, unknown filter '%s'", extra_filters)
        q = Q(pk__isnull=False)  # No op
    return q


def get_node_count_colors(css_property):
    """ Returns a list of tuples with last element being a dict,
        css_property of "color" = [('ClinVar', {color: #ff0000}), etc] """

    node_count_colors = []
    for label, color in BuiltInFilters.COLORS:
        node_count_colors.append((label, {css_property: color}))

    return node_count_colors


def get_node_counts_mine_and_available(analysis):
    node_count_types = analysis.get_node_count_types()

    labels = dict(BuiltInFilters.CHOICES)
    my_choices = [x[0] for x in node_count_types]
    all_choices = [x[0] for x in BuiltInFilters.CHOICES]
    # Needs to stay in order.
    available_choices = []
    for c in all_choices:
        if c not in my_choices:
            available_choices.append(c)

    my_node_counts_list = []
    for node_count in my_choices:
        my_node_counts_list.append({"pk": node_count,
                                    "css_classes": 'node-count-legend-' + node_count,
                                    "description": labels[node_count]})
    available_node_counts_list = []
    for node_count in available_choices:
        available_node_counts_list.append({"pk": node_count,
                                           "css_classes": 'node-count-legend-' + node_count,
                                           "description": labels[node_count]})

    return my_node_counts_list, available_node_counts_list


def get_node_counts_and_labels_dict(node, counts_to_get):

    # Need to do inner query as distinct needs to be applied
    # before aggregate functions
    qs = node.get_queryset(inner_query_distinct=True)
    aggregate_kwargs = {}
    for count_type in counts_to_get:
        if count_type == BuiltInFilters.TOTAL:
            q = None
        else:
            q = get_extra_filters_q(node.analysis.user, node.analysis.genome_build, count_type)
        aggregate_kwargs[count_type] = Count("pk", filter=q)
    return qs.aggregate(**aggregate_kwargs)
