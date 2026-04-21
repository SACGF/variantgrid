from typing import Iterable, Optional, Any, TypeVar, Generic, Type

import sqlparse
from dataclasses_json import DataClassJsonMixin
from django.db import models, connection, transaction
# 970: Added transaction wrapper due to Postgres hanging query
from django.db.models import QuerySet
from django.db.models.enums import TextChoices, IntegerChoices

from library.cache import timed_cache
from library.constants import DAY_SECS


@transaction.atomic
def run_sql(sql, params=None) -> tuple[Any, int]:
    with connection.cursor() as cursor:
        value = cursor.execute(sql, params)  # Remember it only accepts '%s' not %d etc.
        rowcount = cursor.rowcount
        return value, rowcount


def get_postgresql_version() -> str:
    # Few ways to get this, but we'll go with the simpler one:
    # SHOW server_version - '14.12 (Ubuntu 14.12-0ubuntu0.22.04.1)'
    # select version() - 'PostgreSQL 14.12 (Ubuntu 14.12-0ubuntu0.22.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0, 64-bit'
    with connection.cursor() as cursor:
        cursor.execute('SHOW server_version')
        version = cursor.fetchone()[0]
    return version


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
def get_select_from_where_parts_str(sql_str: str) -> tuple[str, str, str]:
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


def get_queryset_select_from_where_parts(qs: QuerySet) -> tuple[str, str, str]:
    """ Returns (select, from, where) """
    sql_str = queryset_to_sql(qs)
    return get_select_from_where_parts_str(sql_str)


def get_queryset_column_names(queryset: QuerySet, extra_columns: list[str]) -> list[str]:
    extra_names = list(queryset.query.extra_select)
    field_names = list(queryset.query.values_select)
    annotation_names = list(queryset.query.annotation_select)  # aggregate_select => annotation_select in Django 1.8

    column_names = extra_names + field_names + annotation_names + extra_columns
    return column_names


def get_cursor_column_names(cursor):
    return [col[0] for col in cursor.description]


def dictfetchall(cursor, column_names: Optional[Iterable[str]] = None) -> list[dict]:
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
        yield from results


def iter_dictfetchall(cursor, column_names: Optional[Iterable[str]] = None) -> dict:
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


T1 = TypeVar('T1', bound=DataClassJsonMixin)


class JSONDataclassField(models.JSONField, Generic[T1]):
    """
    A field type that will return and expect a DataClassJsonMixin type.
    e.g.
    @dataclass
    class Coordinate(DataClassJsonMixin):
        x: int
        y: int
    Will be serialized into the database, and referring to this field will deserialize so you don't have to convert
    from JSONB into your dataclass objects in custom code, it'll just be handled in ORM

    PyCharm isn't smart enough to work out the execpted type, so just add the type hint after declaring using this in a model
    """

    def __init__(self,
                 dataclass_type: Type[T1],
                 illegal_value_result: Any = None,
                 *args,
                 **kwargs):
        """
        :param dataclass_type: The class to convert to/from JSON, should be @dataclass_json or DataClassJsonMixin
        :param illegal_value_result: If we can't convert the JSONB from the database into dataclass_type, return this instead
        """
        self.dataclass_type = dataclass_type
        self.illegal_value_result = illegal_value_result
        super().__init__(*args, **kwargs)

    @property
    def non_db_attrs(self):
        return super().non_db_attrs + ("dataclass_type", "illegal_value_result",)

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        kwargs["dataclass_type"] = self.dataclass_type
        kwargs["illegal_value_result"] = self.illegal_value_result
        return name, path, args, kwargs

    def from_db_value(self, value, expression, connection) -> Optional[T1]:
        if json_obj := super().from_db_value(value, expression, connection):
            try:
                return self.dataclass_type.from_dict(json_obj)
            except ValueError as ve:
                # TODO raise warning
                print(f"Found illegal value in database {ve}")
                return self.illegal_value_result
        return None

    def get_prep_value(self, value: Optional[T1]):
        if hasattr(value, 'to_dict'):
            return value.to_dict()
        return value


TX = TypeVar('TX')


class ChoicesMixin(Generic[TX]):
    """
    While IntegerField, TextField and CharField all allow choices, they'll still return an integer or string.
    TextChoices or IntegerChoices allow you to add additional methods, so using one of the concrete implementations
    below will gice you an instance of the actual Choices instead of the raw string or integer value.

    Lastly you just need to provide the type and if you don't provide choices, they'll be taken from the type

    PyCharm isn't smart enough to work out the execpted type, so just add the type hint after declaring using this in a model
    """

    def __init__(self, choices_type: Type[TX], *args, **kwargs):
        if "choices" not in kwargs:
            kwargs["choices"] = choices_type.choices
        self.choices_type = choices_type
        super().__init__(*args, **kwargs)

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        kwargs["choices_type"] = self.choices_type
        return name, path, args, kwargs

    @property
    def non_db_attrs(self):
        return super().non_db_attrs + ("choices_type",)

    def from_db_value(self, value, expression, connection) -> Optional[TX]:
        if value:
            return self.choices_type(value)
        return None


T2 = TypeVar('T2', bound=TextChoices)


class TextFieldChoices(ChoicesMixin[T2], models.TextField):
    pass


T3 = TypeVar('T3', bound=TextChoices)


class CharFieldChoices(ChoicesMixin[T3], models.CharField):
    pass


T4 = TypeVar('T4', bound=IntegerChoices)


class IntegerFieldChoices(ChoicesMixin[T4], models.IntegerField):
    pass
