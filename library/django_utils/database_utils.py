"""
Raw-SQL helpers: queryset_to_sql and get_queryset_select_from_where_parts turn a QuerySet into SQL
text to embed in COPY / INSERT statements, dictfetchall / iter_db_results read
cursors, sql_delete_qs deletes by a queryset's WHERE without loading rows (dangerous - read it first),
postgres_arrays formats array literals, long_running_sql / get_active_backend_pids / signal_backends /
wait_for_backends_to_stop find and cancel other connections' running queries, pg_settings / get_pg_setting set
and read this connection's run-time settings, and get_table_row_estimates / get_queryset_row_estimate give planner
row counts without scanning.
"""
import contextlib
import json
import logging
import time
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import timedelta
from typing import Any, Optional, TypeVar, Generic, Type, Callable

import sqlparse
from dataclasses_json import DataClassJsonMixin
from django.db import connection, transaction, models

# 970: Added transaction wrapper due to Postgres hanging query
from django.db.models import QuerySet
from django.db.models.lookups import In
from django.db.models.sql.where import NothingNode
from django.db.models.enums import TextChoices, IntegerChoices
from django_json_widget.widgets import JSONEditorWidget

from library.cache import timed_cache
from library.constants import DAY_SECS


@transaction.atomic
def run_sql(sql, params=None) -> tuple[Any, int]:
    with connection.cursor() as cursor:
        value = cursor.execute(sql, params)  # Remember it only accepts '%s' not %d etc.
        rowcount = cursor.rowcount
        return value, rowcount


@dataclass
class RunningQuery:
    pid: int
    duration: timedelta
    query: str
    state: str


def long_running_sql(min_age_in_seconds: int = 30) -> list[RunningQuery]:
    """ Non-idle connections to this database whose current query started more than min_age_in_seconds ago,
        longest first. Idle ones are pooled connections showing their last query, not running anything """
    sql = """
        SELECT pid, now() - query_start AS duration, query, state
        FROM pg_stat_activity
        WHERE datname = current_database() AND state <> 'idle' AND now() - query_start > interval %s
        ORDER BY duration DESC
    """
    with connection.cursor() as cursor:
        cursor.execute(sql, [f"{min_age_in_seconds} seconds"])
        return [RunningQuery(*row) for row in cursor.fetchall()]


def get_active_backend_pids(query_regex: str) -> list[int]:
    """ Other connections to this database currently running a query matching query_regex (case-insensitive) """
    sql = """
        SELECT pid FROM pg_stat_activity
        WHERE datname = current_database() AND pid <> pg_backend_pid() AND state = 'active' AND query ~* %s
    """
    with connection.cursor() as cursor:
        cursor.execute(sql, [query_regex])
        return [row[0] for row in cursor.fetchall()]


def signal_backends(pids: Iterable[int], terminate: bool = False) -> list[int]:
    """ pg_cancel_backend (stop the current query, keep the connection) or, with terminate, pg_terminate_backend
        (close the connection). Returns the pids Postgres signalled - it only signals, so a caller that needs the
        locks released must wait_for_backends_to_stop. Needs the same DB role as the backend (or superuser) """
    pids = list(pids)
    if not pids:
        return []
    pg_function = "pg_terminate_backend" if terminate else "pg_cancel_backend"
    with connection.cursor() as cursor:
        cursor.execute(f"SELECT pid FROM unnest(%s::int[]) AS pid WHERE {pg_function}(pid)", [pids])
        return [row[0] for row in cursor.fetchall()]


def wait_for_backends_to_stop(pids: Iterable[int], timeout_seconds: float) -> bool:
    """ Polls until none of pids is running a query, or timeout_seconds passes. Returns whether they all stopped """
    pids = list(pids)
    if not pids:
        return True

    sql = "SELECT count(*) FROM pg_stat_activity WHERE pid = ANY(%s) AND state = 'active'"
    deadline = time.monotonic() + timeout_seconds
    while True:
        with connection.cursor() as cursor:
            cursor.execute(sql, [pids])
            if not cursor.fetchone()[0]:
                return True
        if time.monotonic() >= deadline:
            logging.warning("Signalled backends (pids=%s) still running after %s seconds", pids, timeout_seconds)
            return False
        time.sleep(0.1)


def get_pg_setting(name: str) -> str:
    """ Current value of a Postgres run-time setting on this connection, as SHOW would print it """
    with connection.cursor() as cursor:
        cursor.execute("SELECT current_setting(%s)", [name])
        return cursor.fetchone()[0]


def _set_pg_settings(values: dict[str, str], local: bool):
    with connection.cursor() as cursor:
        for name, value in values.items():
            cursor.execute("SELECT set_config(%s, %s, %s)", [name, value, local])


@contextlib.contextmanager
def pg_settings(local: bool = False, **values):
    """ Set Postgres run-time settings (eg statement_timeout=5000, work_mem="4GB") on this connection for the
        block. A value of None leaves that setting alone.

        Session settings are put back to what they were on exit - connections are reused (CONN_MAX_AGE) and
        this may be nested inside another caller's settings. local=True is SET LOCAL: it ends with the
        current transaction, so there is nothing to put back (and nothing can run in an aborted one) """
    values = {name: str(value) for name, value in values.items() if value is not None}
    if not values or connection.vendor != 'postgresql':
        yield
        return

    previous = {} if local else {name: get_pg_setting(name) for name in values}
    _set_pg_settings(values, local)
    try:
        yield
    finally:
        if previous:
            _set_pg_settings(previous, local=False)


def get_table_row_estimates(table_names: Optional[Iterable[str]] = None) -> dict[str, int]:
    """ {table: planner row estimate} from pg_class for tables and partitioned tables - all of them, or just
        table_names. -1 means the table has never been vacuumed or analyzed. An inheritance parent only
        counts its own rows, not its children's """
    sql = "SELECT relname, reltuples::bigint FROM pg_class WHERE relkind IN ('r', 'p')"
    params = []
    if table_names is not None:
        sql += " AND relname = ANY(%s)"
        params.append(list(table_names))
    with connection.cursor() as cursor:
        cursor.execute(sql, params)
        return dict(cursor.fetchall())


def get_queryset_row_estimate(qs: QuerySet) -> int:
    """ The planner's row estimate for a queryset, without running it - cheap where count() would scan """
    sql, params = qs.query.sql_with_params()
    with connection.cursor() as cursor:
        cursor.execute(f"EXPLAIN (FORMAT JSON) {sql}", params)
        plan = cursor.fetchone()[0]
    if isinstance(plan, str):
        plan = json.loads(plan)
    return int(plan[0]["Plan"]["Plan Rows"])


def get_postgresql_version() -> str:
    # Few ways to get this, but we'll go with the simpler one:
    # SHOW server_version - '14.12 (Ubuntu 14.12-0ubuntu0.22.04.1)'
    # select version() - 'PostgreSQL 14.12 (Ubuntu 14.12-0ubuntu0.22.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0, 64-bit'
    with connection.cursor() as cursor:
        cursor.execute('SHOW server_version')
        version = cursor.fetchone()[0]
    return version


@contextlib.contextmanager
def render_empty_result_set_sql():
    """ Django short circuits provably empty predicates (eg "pk__in=[]") by raising EmptyResultSet during
        compilation, so there's no SQL to look at. For debugging we want to see the query anyway, so compile
        those to equivalent SQL that matches nothing.

        Patches are global for the duration - a concurrent query that would have short circuited runs the
        (still empty) SQL instead. """

    def _in_process_rhs(self, compiler, connection_):
        if self.rhs_is_direct_value() and not [r for r in self.rhs if r is not None]:
            return "(NULL)", []
        return original_in_process_rhs(self, compiler, connection_)

    def _nothing_as_sql(*_args, **_kwargs):
        return "0 = 1", []

    original_in_process_rhs = In.process_rhs
    original_nothing_as_sql = NothingNode.as_sql
    In.process_rhs = _in_process_rhs
    NothingNode.as_sql = _nothing_as_sql
    try:
        yield
    finally:
        In.process_rhs = original_in_process_rhs
        NothingNode.as_sql = original_nothing_as_sql


def queryset_to_sql(queryset: QuerySet, pretty=False) -> str:
    """ str(queryset.query) doesn't quote variables properly....

        From: https://stackoverflow.com/a/47542953
        qs.query returns something that isn't valid SQL, this returns the actual
        valid SQL that's executed: https://code.djangoproject.com/ticket/17741  """

    query, params = queryset.query.sql_with_params()
    PREFIX = 'select 1 -- '
    with connection.cursor() as cursor:
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


def sql_delete_qs(qs, batch_size: Optional[int] = None) -> int:
    """ ** WARNING DANGEROUS ***
        A way to perform deletes in batches, in the DB
        returns rows deleted """

    pk_qs = qs.values_list("pk", flat=True)
    meta = qs.model._meta
    if batch_size:
        limit = f"LIMIT {int(batch_size)}"
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
    DO NOT USE - doesn't work well with auditlog which double encodes, likely a problem with serialize as well
    """

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

    def get_fallback_value(self):
        fallback_value = self.illegal_value_result
        if fallback_value is None:
            fallback_value = self.default
        if isinstance(fallback_value, Callable):
            fallback_value = fallback_value()
        if isinstance(fallback_value, dict):
            fallback_value = self.dataclass_type.from_dict(fallback_value)
        return fallback_value

    # def from_db_value(self, value, expression, connection) -> Optional[T1]:
    #     if json_obj := super().from_db_value(value, expression, connection):
    #         try:
    #             return self.dataclass_type.from_dict(json_obj)
    #         except Exception as ex:
    #             # TODO raise warning
    #             print(f"Found illegal value in database {value.__class__} \"{value}\"")
    #             self.get_fallback_value()
    #     return None

    def from_db_value(self, value, expression, connection):
        """Converts database JSON string/object to a dataclass instance."""
        if value is None:
            return None
        try:
            if isinstance(value, dict):
                return self.dataclass_type.from_dict(value)
            if isinstance(value, str):
                return self.dataclass_type.from_json(value)
        except:
            print(f"Could not convert value from DB {value.__class__} \"{value}\" to {self.dataclass_type}")
        return value

    def get_prep_value(self, value: Optional[T1]):
        if value is None:
            return None
        if hasattr(value, 'to_dict'):
            return value.to_dict()
        if not isinstance(value, dict):
            # print(f"WARNING - prep value doesn't have to_dict or isn't a dict, it's {value.__class__}")
            raise ValueError(f"WARNING - prep value doesn't have to_dict or isn't a dict, it's {value.__class__}")
        return value


TX = TypeVar('TX')


class ChoicesMixin(Generic[TX]):
    """
    While IntegerField, TextField and CharField all allow choices, they'll still return an integer or string.
    TextChoices or IntegerChoices allow you to add additional methods, so using one of the concrete implementations
    below will give you an instance of the actual Choices instead of the raw string or integer value.

    Lastly you just need to provide the type and if you don't provide choices, they'll be taken from the type

    PyCharm isn't smart enough to work out the expected type, so just add the type hint after declaring using this in a model
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
        if value is not None:
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


class JSONDataClassAdminWidget(JSONEditorWidget):

    def format_value(self, value):
        if not isinstance(value, (dict, list)):
            value = json.loads(value)
            if isinstance(value, str):
                # TODO, why would these values be double escaped anyway?
                value = json.loads(value)
            return value
        else:
            return value
