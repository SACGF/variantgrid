from django.db import connections
from django.db.models.query import QuerySet
from django.db.models.sql import UpdateQuery
from django.db.models.sql.compiler import SQLCompiler, SQLUpdateCompiler
from django.db.models.sql.query import Query


def get_queryset_with_transformer_hook(klass):
    sql_transformers = []

    # Have to make a whole bunch of classes just so that the compiler can be passed in transformers
    # Table with the appropriate partition
    class TransformerSQLMixin:

        def as_sql(self, **kwargs):
            sql, params = super().as_sql(**kwargs)
            for transformer in sql_transformers:
                sql = transformer(sql)
            return sql, params

    class TransformerSQLCompiler(TransformerSQLMixin, SQLCompiler):
        pass

    class TransformerSQLUpdateCompiler(TransformerSQLMixin, SQLUpdateCompiler):
        pass

    class QueryTransformerCompilerMixin:
        compiler_klass = None

        def get_compiler(self, using=None, connection=None):
            if using is None and connection is None:
                raise ValueError("Need either using or connection")
            if using:
                connection = connections[using]
            return self.compiler_klass(self, connection, using)

    class TransformerUpdateQuery(QueryTransformerCompilerMixin, UpdateQuery):
        compiler_klass = TransformerSQLUpdateCompiler

    class TransformerQuery(QueryTransformerCompilerMixin, Query):
        compiler_klass = TransformerSQLCompiler

        def chain(self, klass=None):
            if klass == UpdateQuery:
                klass = TransformerUpdateQuery
            return super().chain(klass=klass)

    class TransformerQuerySet(QuerySet):

        def __init__(self, model=None, query=None, using=None, hints=None):
            # print(f"Queryset({model}, {query}, {using}, {hints})")
            if query is None:
                query = TransformerQuery(model)
            super().__init__(model, query, using, hints)

        def add_sql_transformer(self, transformer):
            sql_transformers.append(transformer)

    qs = TransformerQuerySet(model=klass)
    return qs.all()
