"""
Helper functions for grid to be able to use custom SQL to join to sample columns

This is in snpdb not analysis as it needs to be called by snpdb.tasks.cohort_genotype_tasks
and I don't want cross-dependencies between those apps
"""
from collections import defaultdict

from library.utils.database_utils import get_queryset_select_from_where_parts
from library.utils import single_quote
from snpdb.models import CohortGenotype, VCF


def get_left_outer_join_on_variant(partition_table):
    return f'LEFT OUTER JOIN "{partition_table}" ON ("snpdb_variant"."id" = "{partition_table}"."variant_id")'


def get_available_format_columns(cohorts):
    available_format_columns = {
        # We want to always show some fields
        "samples_zygosity": True,
        "samples_allele_depth": True,
        "samples_allele_frequency": True,
        "samples_read_depth": True,
        "samples_genotype_quality": False,
        "samples_phred_likelihood": False,
        "samples_filters": False,
    }

    vcf_qs = VCF.objects.filter(sample__cohortsample__cohort__in=cohorts).distinct()
    for gq, pl, ft in vcf_qs.values_list("genotype_quality_field", "phred_likelihood_field", "sample_filters_field"):
        if gq:
            available_format_columns["samples_genotype_quality"] = True
        if pl:
            available_format_columns["samples_phred_likelihood"] = True
        if ft:
            available_format_columns["samples_filters"] = True
    return available_format_columns


def get_columns_and_sql_parts_for_cohorts(values_queryset, cohorts):
    new_column_pack_lists = defaultdict(list)
    new_joins = []
    available_format_columns = get_available_format_columns(cohorts)

    for cohort in cohorts:
        partition_table = cohort.cohort_genotype_collection.get_partition_table()
        for column, (is_array, empty_value) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
            if not available_format_columns[column]:
                continue  # No data, so don't show

            if is_array:
                array_values = ", ".join(str(s) for s in [empty_value] * cohort.sample_count)
                empty_data = f"array[{array_values}]"
            else:
                empty_data = empty_value * cohort.sample_count
                empty_data = single_quote(empty_data)

            pack_entry = f"coalesce({partition_table}.{column}, {empty_data})"
            new_column_pack_lists[column].append(pack_entry)
        new_joins.append(get_left_outer_join_on_variant(partition_table))

    return _get_columns_and_sql_parts_for_pack_lists_and_joins(values_queryset, new_column_pack_lists, new_joins)


def _get_columns_and_sql_parts_for_pack_lists_and_joins(values_queryset, new_column_pack_lists, new_joins):
    select_part, from_part, where_part = get_queryset_select_from_where_parts(values_queryset)
    new_columns = []
    new_select = []
    for column, pack_entry in new_column_pack_lists.items():
        packed_column = f"packed_{column}"
        packed_source = " || ".join(pack_entry)  # || works for both string and arrays
        new_select.append(f"{packed_source} as {packed_column}")
        new_columns.append(packed_column)

    select_part = ",\n".join([select_part] + new_select)
    from_part = "\n".join([from_part] + new_joins)
    return new_columns, select_part, from_part, where_part
