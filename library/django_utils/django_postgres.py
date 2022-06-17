from typing import List

from django.db import models
from django.db.models import Model, sql


class PostgresRealField(models.Field):
    """ 32 bit float """
    def db_type(self, connection):
        return 'real'


def pg_sql_array(values):
    return 'array[%s]' % ','.join(values)


def _escape_sql_param(param):
    if param is None:
        param = 'NULL'
    elif isinstance(param, List):
        param = pg_sql_array(map(str, param))
    elif isinstance(param, str):
        return f"'{param}'"
    return str(param)


# noinspection PyProtectedMember
def model_to_insert_sql(model_list: List[Model], db_table: str = None, ignore_fields: List[str] = None):
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
    setattr(compiler, 'return_id', False)
    raw_statements = compiler.as_sql()
    model._meta.db_table = old_table  # Put table name back

    mixed_statements = [statement % tuple(_escape_sql_param(param) for param in params)
                        for statement, params in raw_statements]
    return mixed_statements
