import csv
import os
import shutil
import tempfile

from django.db import connection
from django.test import TestCase

from library.django_utils.database_utils import postgres_array_quote, postgres_arrays
from library.utils.file_utils import IteratorFile
from upload.vcf.sql_copy_files import sql_copy_csv, sql_copy_csv_file, write_sql_copy_csv

TEMP_TABLE = "sql_copy_test"
COLUMNS = ["a", "b"]


class SqlCopyFilesTest(TestCase):

    def setUp(self):
        with connection.cursor() as cursor:
            cursor.execute(f'CREATE TEMP TABLE {TEMP_TABLE} (a int, b text, "end" int) ON COMMIT DROP')

    def _get_rows(self):
        with connection.cursor() as cursor:
            cursor.execute(f"SELECT a, b FROM {TEMP_TABLE} ORDER BY a")
            return cursor.fetchall()

    def _write_csv(self, contents) -> str:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(contents)
            self.addCleanup(os.unlink, f.name)
            return f.name

    def test_sql_copy_csv(self):
        filename = self._write_csv("1,one\n2,\n")
        self.assertEqual(sql_copy_csv(filename, TEMP_TABLE, COLUMNS), 2)
        self.assertEqual(self._get_rows(), [(1, "one"), (2, None)])

    def test_sql_copy_csv_tab_delimiter(self):
        filename = self._write_csv("1\tone\n")
        self.assertEqual(sql_copy_csv(filename, TEMP_TABLE, COLUMNS, delimiter='\t'), 1)
        self.assertEqual(self._get_rows(), [(1, "one")])

    def test_sql_copy_csv_quoted(self):
        filename = self._write_csv('1,"one,two"\n')
        self.assertEqual(sql_copy_csv(filename, TEMP_TABLE, COLUMNS, quote='"'), 1)
        self.assertEqual(self._get_rows(), [(1, "one,two")])

    def test_sql_copy_csv_reserved_word_column(self):
        """ VARIANTS_HEADER has 'end' in it, which needs quoting """

        filename = self._write_csv("1,one,42\n")
        self.assertEqual(sql_copy_csv(filename, TEMP_TABLE, COLUMNS + ["end"]), 1)
        with connection.cursor() as cursor:
            cursor.execute(f'SELECT "end" FROM {TEMP_TABLE}')
            self.assertEqual(cursor.fetchall(), [(42,)])

    def test_text_array_round_trip(self):
        """ CohortGenotype.samples_filters: array elements quoted as BulkGenotypeVCFProcessor does, written
            with its CSV options - the '\\' filter code must survive both COPY and array parsing #2480 """

        values = ["\\", "\\a", '"', "a,b", "{}"]
        array_str = postgres_arrays([*(postgres_array_quote(v) for v in values), None])
        filename = os.path.join(tempfile.mkdtemp(), "array.csv")
        self.addCleanup(shutil.rmtree, os.path.dirname(filename))
        write_sql_copy_csv([[1, array_str]], filename, quoting=csv.QUOTE_NONE, escapechar='\\')

        with connection.cursor() as cursor:
            cursor.execute("CREATE TEMP TABLE sql_copy_array_test (a int, b text[]) ON COMMIT DROP")
            sql_copy_csv(filename, "sql_copy_array_test", COLUMNS)
            cursor.execute("SELECT b FROM sql_copy_array_test")
            self.assertEqual(cursor.fetchone()[0], [*values, None])

    def test_sql_copy_csv_file_iterator(self):
        """ BulkVCFCountInserter streams rows in via IteratorFile rather than a file on disk """

        f = IteratorFile(iter(["1,one\n", "2,two\n"]))
        self.assertEqual(sql_copy_csv_file(f, TEMP_TABLE, COLUMNS), 2)
        self.assertEqual(self._get_rows(), [(1, "one"), (2, "two")])
