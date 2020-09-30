"""
Helper functions for grid to be able to use custom SQL to join to sample columns

This is in snpdb not analysis as it needs to be called by snpdb.tasks.cohort_genotype_tasks
and I don't want cross-dependencies between those apps
"""
from collections import defaultdict

from library.database_utils import get_queryset_select_from_where_parts
from library.utils import single_quote
from snpdb.models import CohortGenotype


def get_left_outer_join_on_variant(partition_table):
    return f'LEFT OUTER JOIN "{partition_table}" ON ("snpdb_variant"."id" = "{partition_table}"."variant_id")'


def get_columns_and_sql_parts_for_cohorts(values_queryset, cohorts):
    new_column_pack_lists = defaultdict(list)
    new_joins = []

    for cohort in cohorts:
        num_samples = cohort.get_samples().count()
        partition_table = cohort.cohort_genotype_collection.get_partition_table()
        for column, (is_array, empty_value) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
            if is_array:
                array_values = ", ".join(str(s) for s in [empty_value] * num_samples)
                empty_data = f"array[{array_values}]"
            else:
                empty_data = empty_value * num_samples
                empty_data = single_quote(empty_data)

            pack_entry = f"coalesce({partition_table}.{column}, {empty_data})"
            new_column_pack_lists[column].append(pack_entry)
        new_joins.append(get_left_outer_join_on_variant(partition_table))

    return get_columns_and_sql_parts_for_pack_lists_and_joins(values_queryset, new_column_pack_lists, new_joins)


def get_columns_and_sql_parts_for_pack_lists_and_joins(values_queryset, new_column_pack_lists, new_joins):
    select_part, from_part, where_part = get_queryset_select_from_where_parts(values_queryset)
    new_columns = []
    new_select = []
    for column in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE:
        packed_column = f"packed_{column}"
        packed_source = " || ".join(new_column_pack_lists[column])  # || works for both string and arrays
        new_select.append(f"{packed_source} as {packed_column}")
        new_columns.append(packed_column)

    select_part = ",\n".join([select_part] + new_select)

    # The SQL spec states that explicit joins are performed before implicit joins.
    # Django will already manage this for us, but we have to keep it up, so insert our new joins at the start
    from_before_comma, *existing_joins = from_part.split(',')
    if existing_joins:
        print("*" * 50)
        print("existing_joins:")
        print(existing_joins)
    from_before_comma = "\n".join([from_before_comma] + new_joins)
    from_part = ','.join([from_before_comma] + existing_joins)
    return new_columns, select_part, from_part, where_part
