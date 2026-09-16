
from django.db import connection, models
from django.db.backends.postgresql.psycopg_any import mogrify
from django.db.models import Model, sql

# IteratorFile.read() concatenates until it has this much, so big reads go quadratic
COPY_READ_SIZE = 8192


class PostgresRealField(models.Field):
    """ 32 bit float """
    def db_type(self, connection):
        return 'real'


# psycopg2 and psycopg3 differ on the driver APIs below, and deployments can be on either, so go
# through these wrappers rather than the raw connection/cursor. Both drivers are tested in
# library.django_utils.tests.test_django_postgres


def get_backend_pid(cursor) -> int:
    """ PID of the Postgres backend serving this cursor's connection.
        psycopg3 dropped get_backend_pid() - connection.info works on both drivers """
    return cursor.db.connection.info.backend_pid


def copy_from_file(cursor, copy_sql: str, f):
    """ Runs a 'COPY <table> FROM STDIN' statement, streaming from a file-like object.
        psycopg3 replaced copy_expert()/copy_from() with the cursor.copy() context manager """
    if hasattr(cursor, "copy_expert"):  # psycopg2
        cursor.copy_expert(copy_sql, f)
    else:
        with cursor.copy(copy_sql) as copy:
            while data := f.read(COPY_READ_SIZE):
                copy.write(data)
    return cursor.rowcount


def pg_sql_array(values):
    return f'array[{",".join(values)}]'


# noinspection PyProtectedMember
def model_to_insert_sql(model_list: list[Model], db_table: str = None, ignore_fields: list[str] = None):
    """ From https://stackoverflow.com/a/63715608/295724 """

    if ignore_fields is None:
        ignore_fields = []

    model = model_list[0]
    old_table = model._meta.db_table
    if db_table:
        model._meta.db_table = db_table
    fields = [f for f in model._meta.local_fields if f.name not in ignore_fields]
    q = sql.InsertQuery(model)
    q.insert_values(fields, model_list)

    compiler = q.get_compiler('default')
    # Normally, execute sets this, but we don't want to call execute
    compiler.return_id = False
    raw_statements = compiler.as_sql()
    model._meta.db_table = old_table  # Put table name back

    # mogrify does the client side param binding - it handles the driver specific wrappers Django puts
    # around params (eg Jsonb for JSONField) as well as quoting
    return [mogrify(statement, params, connection) for statement, params in raw_statements]
