from django.db import connection, transaction
import sqlparse


#970: Added transaction wrapper due to Postgres hanging query
@transaction.atomic
def run_sql(sql, params=None):
    with connection.cursor() as cursor:
        value = cursor.execute(sql, params)  # Remember it only accepts '%s' not %d etc.
        rowcount = cursor.rowcount
        return value, rowcount


def queryset_to_sql(queryset, pretty=False):
    """ str(queryset.query) doesn't quote variables properly....

        From: https://stackoverflow.com/a/47542953
        qs.query returns something that isn't valid SQL, this returns the actual
        valid SQL that's executed: https://code.djangoproject.com/ticket/17741  """

    cursor = connection.cursor()
    query, params = queryset.query.sql_with_params()
    cursor.execute('EXPLAIN ' + query, params)
    res = str(cursor.db.ops.last_executed_query(cursor, query, params))
    assert res.startswith('EXPLAIN ')
    query_sql = res[len('EXPLAIN '):]

    if pretty:
        query_sql = sqlparse.format(query_sql, reindent=True, keyword_case='upper')

    return query_sql


def get_select_from_where_parts_str(sql_str):
    from_pos = sql_str.find("FROM")
    where_pos = sql_str.find("WHERE")
    if where_pos < 0:
        where_pos = len(sql_str)

    select_part = sql_str[:from_pos]
    from_part = sql_str[from_pos:where_pos]
    where_part = sql_str[where_pos:]
    return select_part, from_part, where_part


def get_queryset_select_from_where_parts(qs):
    """ Returns (select, from, where) """
    sql_str = queryset_to_sql(qs)
    return get_select_from_where_parts_str(sql_str)


def get_queryset_column_names(queryset, extra_columns):
    extra_names = list(queryset.query.extra_select)
    field_names = list(queryset.query.values_select)
    annotation_names = list(queryset.query.annotation_select)  # aggregate_select => annotation_select in Django 1.8

    column_names = extra_names + field_names + annotation_names + extra_columns
    return column_names


def get_cursor_column_names(cursor):
    return [col[0] for col in cursor.description]


def dictfetchall(cursor, column_names=None):
    if column_names is None:
        column_names = get_cursor_column_names(cursor)

    return [dict(list(zip(column_names, row))) for row in cursor.fetchall()]


# From http://code.activestate.com/recipes/137270-use-generators-for-fetching-large-db-record-sets/
def iter_db_results(cursor, arraysize=1000):
    'An iterator that uses fetchmany to keep memory usage down'
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for tup in results:
            yield tup


def iter_dictfetchall(cursor, column_names=None):
    if column_names is None:
        column_names = get_cursor_column_names(cursor)

    for row in iter_db_results(cursor, 10000):
        data = dict(list(zip(column_names, row)))
        yield data


def sql_delete_qs(qs, batch_size=None):
    """ ** WARNING DANGEROUS ***
        A way to perform deletes in batches, in the DB
        returns rows deleted """

    pk_qs = qs.values_list("pk", flat=True)
    meta = qs.model._meta
    if batch_size:
        limit = f"LIMIT {batch_size}"
    else:
        limit = ""

    sql = f"DELETE FROM {meta.db_table} where {meta.pk.name} in ({pk_qs.query} {limit})"
    total_rowcount = 0
    rowcount = True
    while rowcount:
        _, rowcount = run_sql(sql)
        total_rowcount += rowcount
    return total_rowcount
