from typing import Iterable, Optional, Tuple, List, Any, Dict

import sqlparse
from django.db import connection, transaction
# 970: Added transaction wrapper due to Postgres hanging query
from django.db.models import QuerySet

from library.cache import timed_cache
from library.constants import DAY_SECS


@transaction.atomic
def run_sql(sql, params=None) -> Tuple[Any, int]:
    with connection.cursor() as cursor:
        value = cursor.execute(sql, params)  # Remember it only accepts '%s' not %d etc.
        rowcount = cursor.rowcount
        return value, rowcount


def queryset_to_sql(queryset: QuerySet, pretty=False) -> str:
    """ str(queryset.query) doesn't quote variables properly....

        From: https://stackoverflow.com/a/47542953
        qs.query returns something that isn't valid SQL, this returns the actual
        valid SQL that's executed: https://code.djangoproject.com/ticket/17741  """

    cursor = connection.cursor()
    query, params = queryset.query.sql_with_params()
    PREFIX = 'select 1 -- '
    cursor.execute(PREFIX + query, params)
    res = str(cursor.db.ops.last_executed_query(cursor, query, params))
    assert res.startswith(PREFIX)
    query_sql = res[len(PREFIX):]

    if pretty:
        query_sql = sqlparse.format(query_sql, reindent=True, keyword_case='upper')

    if queryset.query.has_select_fields:
        assert query_sql.upper().startswith("SELECT"), "Select query startswith SELECT"

    return query_sql


@timed_cache(ttl=DAY_SECS)
def get_select_from_where_parts_str(sql_str: str) -> Tuple[str, str, str]:
    parsed = sqlparse.parse(sql_str)
    tokens = parsed[0].tokens
    from_token_index = None
    where_token_index = None

    for i, token in enumerate(tokens):
        if token.is_keyword:
            if token.value.upper() == "FROM":
                from_token_index = i
        elif isinstance(token, sqlparse.sql.Where):
            where_token_index = i

    if where_token_index is None:
        where_token_index = len(tokens)

    select_statement = sqlparse.sql.Statement(tokens[:from_token_index])
    from_statement = sqlparse.sql.Statement(tokens[from_token_index:where_token_index])
    where_statement = sqlparse.sql.Statement(tokens[where_token_index:])
    return str(select_statement), str(from_statement), str(where_statement)


def get_queryset_select_from_where_parts(qs: QuerySet) -> Tuple[str, str, str]:
    """ Returns (select, from, where) """
    sql_str = queryset_to_sql(qs)
    return get_select_from_where_parts_str(sql_str)


def get_queryset_column_names(queryset: QuerySet, extra_columns: List[str]) -> List[str]:
    extra_names = list(queryset.query.extra_select)
    field_names = list(queryset.query.values_select)
    annotation_names = list(queryset.query.annotation_select)  # aggregate_select => annotation_select in Django 1.8

    column_names = extra_names + field_names + annotation_names + extra_columns
    return column_names


def get_cursor_column_names(cursor):
    return [col[0] for col in cursor.description]


def dictfetchall(cursor, column_names: Optional[Iterable[str]] = None) -> List[Dict]:
    if column_names is None:
        column_names = get_cursor_column_names(cursor)

    return [dict(list(zip(column_names, row))) for row in cursor.fetchall()]


# From http://code.activestate.com/recipes/137270-use-generators-for-fetching-large-db-record-sets/
def iter_db_results(cursor, array_size=1000):
    """ An iterator that uses fetchmany to keep memory usage down """
    while True:
        results = cursor.fetchmany(array_size)
        if not results:
            break
        for tup in results:
            yield tup


def iter_dictfetchall(cursor, column_names: Optional[Iterable[str]] = None) -> Dict:
    if column_names is None:
        column_names = get_cursor_column_names(cursor)

    for row in iter_db_results(cursor, 10000):
        yield dict(zip(column_names, row))


def sql_delete_qs(qs, batch_size: Optional[int] = None) -> int:
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


def postgres_arrays(array):
    return "{%s}" % ','.join([str(s) if s is not None else "NULL" for s in array])
