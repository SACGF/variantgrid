import os
import tempfile

from django.db import connection
from django.test import TestCase

from library.utils.file_utils import IteratorFile
from upload.vcf.sql_copy_files import sql_copy_csv, sql_copy_csv_file

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

    def test_sql_copy_csv_file_iterator(self):
        """ BulkVCFCountInserter streams rows in via IteratorFile rather than a file on disk """

        f = IteratorFile(iter(["1,one\n", "2,two\n"]))
        self.assertEqual(sql_copy_csv_file(f, TEMP_TABLE, COLUMNS), 2)
        self.assertEqual(self._get_rows(), [(1, "one"), (2, "two")])
