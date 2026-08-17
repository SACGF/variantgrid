from django.db import connection
from django.test import TestCase

from library.django_utils.django_postgres import copy_from_file, get_backend_pid
from library.utils.file_utils import IteratorFile


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
