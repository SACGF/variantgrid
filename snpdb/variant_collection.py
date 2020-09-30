from django.db import connection
import logging

# TODO: Could possibly use QS ops - https://docs.djangoproject.com/en/2.0/ref/models/querysets/#union
# Rewritten due to issue #221 - Original implementation had problems when
# Joining to a table with multiple entries for a variant (the OR caused a cartesian product)


def write_sql_to_variant_collection(variant_collection, sql):
    table_name = variant_collection.get_partition_table()

    sql_template = """insert into %(table_name)s
                     (variant_id, variant_collection_id)
                     %(sql)s;"""
    sql = sql_template % {'table_name': table_name, 'sql': sql}
    # print(sql)

    cursor = connection.cursor()
    cursor.execute(sql)
    variant_collection.count = cursor.rowcount
    variant_collection.save()
    logging.debug("write_sql_to_variant_collection: count = %d", variant_collection.count)


def write_variant_set_operation(variant_collection, a_sql, b_sql, set_operation):
    # Use SQL ops as (a.get_queryset() | b.get_queryset()).distinct() is slow
    sql_template = """select id, %(variant_collection_id)d FROM (
             %(a_sql)s
             %(set_operation)s
             %(b_sql)s
             ) I"""
    sql = sql_template % {'a_sql': a_sql,
                          'b_sql': b_sql,
                          'variant_collection_id': variant_collection.pk,
                          'set_operation': set_operation}
    write_sql_to_variant_collection(variant_collection, sql)
