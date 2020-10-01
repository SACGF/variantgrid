from collections import defaultdict

from django.contrib.auth.models import User
from django.core.exceptions import EmptyResultSet
from django.db import connection
from django.db.models.query_utils import Q
import logging

from annotation.models.damage_enums import PathogenicityImpact
from library.database_utils import get_queryset_select_from_where_parts, dictfetchall
from snpdb.models.models_enums import BuiltInFilters
from classification.enums import ClinicalSignificance
from classification.models import Classification, GenomeBuild

# Add the necessary fields to qs to create join:
REQUIRED_FIELDS = [
    "clinvar__highest_pathogenicity",
    "variantannotation__gene__ensemblgeneannotation__omim_phenotypes",
    "variantannotation__impact"
]

CLASSIFICATION_COUNT_SQL = """
select 1
from classification_classification
where classification_classification.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where allele_id in (
        select allele_id
        from snpdb_variantallele
        where variant_id = snpdb_variant.id
    )
)
"""

COUNTS = {
    BuiltInFilters.TOTAL: "count(*)",
    BuiltInFilters.CLINVAR: "sum(case when %(annotation_clinvar)s.highest_pathogenicity >= 4 then 1 else 0 end)",
    BuiltInFilters.OMIM: "sum(case when %(annotation_ensemblgeneannotation)s.omim_phenotypes is not null then 1 else 0 end)",
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: "sum(case when %(annotation_variantannotation)s.impact in ('H', 'M') then 1 else 0 end)",
    BuiltInFilters.COSMIC: "sum(case when %(annotation_variantannotation)s.cosmic_id is not null then 1 else 0 end)",
    BuiltInFilters.CLASSIFIED: f"sum(case when exists ({CLASSIFICATION_COUNT_SQL}) then 1 else 0 end)",
    BuiltInFilters.CLASSIFIED_PATHOGENIC: f"sum(case when exists ({CLASSIFICATION_COUNT_SQL} AND classification_classification.clinical_significance in ('4', '5')) then 1 else 0 end)"
}


SELECT_INTERNALLY_CLASSIFIED_SQL = """
select string_agg(coalesce(classification_classification.clinical_significance, 'U'), '|')
from classification_classification
where classification_classification.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where
    allele_id in (
        select allele_id from snpdb_variantallele where variant_id = snpdb_variant.id
    )
)
"""

SELECT_MAX_INTERNAL_CLASSIFICATION = """
select max(coalesce(classification_classification.clinical_significance, '0'))
from classification_classification
where classification_classification.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where
    allele_id in (
        select allele_id from snpdb_variantallele where variant_id = snpdb_variant.id
    )
)
"""

INTERNAL_CLASSIFICATION_ALIASES_AND_SELECT = {
    "internally_classified": SELECT_INTERNALLY_CLASSIFIED_SQL,
    "max_internal_classification": SELECT_MAX_INTERNAL_CLASSIFICATION,
}


def get_extra_filters_q(user: User, genome_build: GenomeBuild, extra_filters):
    if extra_filters == BuiltInFilters.CLINVAR:
        q = Q(clinvar__highest_pathogenicity__gte=4)
    elif extra_filters == BuiltInFilters.OMIM:
        q = Q(variantannotation__gene__ensemblgeneannotation__omim_phenotypes__isnull=False)
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
    for (label, color) in BuiltInFilters.COLORS:
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


def get_node_counts_and_labels_dict(node):

    # TODO: We should pass in the labels we want, only join to the appropriate tables and retrieve what we want
    # so if we only want clinvar or classified we only have to scan short tables

    # Need to do inner query as distinct needs to be applied
    # before aggregate functions
    qs = node.get_queryset(inner_query_distinct=True)
    qs = qs.values(*REQUIRED_FIELDS)

    def get_count_alias(count_type):
        return f"{count_type}_count".lower()

    try:
        _, from_str, where_str = get_queryset_select_from_where_parts(qs)

        partition_names = node.analysis.annotation_version.get_partition_names()

        select_columns = []
        for (count_type, column_string) in COUNTS.items():
            column_string %= partition_names
            column_string += " as " + get_count_alias(count_type)
            select_columns.append(column_string)

        select_str = 'SELECT ' + ',\n'.join(select_columns)
        sql = '\n'.join([select_str, from_str, where_str])
        # logging.info("NODE COUNT sql was:")
        # logging.info(sql)

        try:
            cursor = connection.cursor()
            cursor.execute(sql)
        except Exception as e:
            logging.error(e)
            logging.error(sql)
            raise

        data = dictfetchall(cursor)
        if len(data) != 1:
            msg = f"Expected single row! Was {len(data)} rows"
            raise ValueError(msg)
        data = data[0]
    except EmptyResultSet:
        data = defaultdict(int)

    node_counts = {}
    for count_type in COUNTS:
        count_alias = get_count_alias(count_type)
        node_counts[count_type] = data[count_alias] or 0

    return node_counts
