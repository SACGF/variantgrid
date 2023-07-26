import logging

from django.db import connection


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
