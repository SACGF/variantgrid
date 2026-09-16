from django.db import connection
from django.test import TestCase

from library.django_utils.django_postgres import copy_from_file, get_backend_pid, model_to_insert_sql
from library.utils.file_utils import IteratorFile
from library.utils.database_utils import run_sql
from snpdb.models.models_clingen_allele import ClinGenAllele


class PsycopgDriverTest(TestCase):
    """ These wrap the psycopg APIs that changed between psycopg2 and psycopg3, so run them
        against whatever driver this deployment has.

        To check the psycopg2 path on a psycopg3 box, put a psycopg.py that raises ImportError
        on the PYTHONPATH - Django then falls back to psycopg2 """

    def _create_temp_table(self, cursor):
        cursor.execute("CREATE TEMP TABLE copy_test (a int, b text) ON COMMIT DROP")

    def test_get_backend_pid(self):
        with connection.cursor() as cursor:
            cursor.execute("SELECT pg_backend_pid()")
            expected_pid = cursor.fetchone()[0]
            self.assertEqual(get_backend_pid(cursor), expected_pid)

    def test_copy_from_file(self):
        """ Text format - empty fields are NULL """

        with connection.cursor() as cursor:
            self._create_temp_table(cursor)
            sql = "COPY copy_test (a,b) FROM STDIN WITH (FORMAT text, DELIMITER ',', NULL '')"
            f = IteratorFile(iter(["1,one\n", "2,\n"]))
            self.assertEqual(copy_from_file(cursor, sql, f), 2, "COPY rowcount")

            cursor.execute("SELECT a, b FROM copy_test ORDER BY a")
            self.assertEqual(cursor.fetchall(), [(1, "one"), (2, None)])

    def test_copy_from_file_csv_quoted(self):
        """ CSV format - quoted fields can contain the delimiter """

        with connection.cursor() as cursor:
            self._create_temp_table(cursor)
            sql = 'COPY copy_test (a,b) FROM STDIN WITH (FORMAT csv, QUOTE \'"\')'
            f = IteratorFile(iter(['1,"one,two"\n']))
            self.assertEqual(copy_from_file(cursor, sql, f), 1, "COPY rowcount")

            cursor.execute("SELECT a, b FROM copy_test")
            self.assertEqual(cursor.fetchall(), [(1, "one,two")])


class ModelToInsertSqlTest(TestCase):
    def test_insert_json_and_quoted_text(self):
        """ JSONField params arrive as driver specific Jsonb wrappers, and text can contain quotes """

        api_response = {"AD": [1, 2], "name": "O'Brien \"quoted\"", "nested": {"pct": 1.5, "ok": True}}
        cga = ClinGenAllele(pk=123456789, api_response=api_response)
        for insert_sql in model_to_insert_sql([cga]):
            run_sql(insert_sql)

        from_db = ClinGenAllele.objects.get(pk=123456789)
        self.assertEqual(from_db.api_response, api_response)

    def test_insert_into_alternate_table(self):
        with connection.cursor() as cursor:
            cursor.execute("CREATE TEMP TABLE clingen_copy (LIKE snpdb_clingenallele) ON COMMIT DROP")

        cga = ClinGenAllele(pk=987654321, api_response={"errorType": "ServerError"})
        for insert_sql in model_to_insert_sql([cga], db_table="clingen_copy"):
            run_sql(insert_sql)

        self.assertFalse(ClinGenAllele.objects.filter(pk=987654321).exists(), "Inserted into other table")
        with connection.cursor() as cursor:
            cursor.execute("SELECT id FROM clingen_copy")
            self.assertEqual(cursor.fetchall(), [(987654321,)])
