"""
Tests for RelatedModelsPartitionModel partition creation.
"""

from django.db import connection
from django.test import TestCase

from genes.models import GeneCoverageCollection
from snpdb.models import DataState, GenomeBuild


class PartitionIdSequenceTests(TestCase):
    """ Children of an inheritance parent get neither the parent's identity property nor a default,
        so the partition has to be pointed at the parent's sequence itself. """

    BASE_TABLE_NAME = "genes_genecoverage"

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch38()

    def _make_gcc(self) -> GeneCoverageCollection:
        return GeneCoverageCollection.objects.create(path="/tmp/nonexistent.tsv",
                                                     data_state=DataState.COMPLETE,
                                                     genome_build=self.genome_build)

    def _assert_partition_id_auto_increments(self, gcc: GeneCoverageCollection):
        """ Evaluate the partition's id default - COPY doesn't provide ids, so it has to hand out its own """

        partition_table = gcc.get_partition_table(base_table_name=self.BASE_TABLE_NAME)
        with connection.cursor() as cursor:
            cursor.execute("""
                SELECT pg_get_expr(d.adbin, d.adrelid)
                FROM pg_attribute a
                JOIN pg_attrdef d ON d.adrelid = a.attrelid AND d.adnum = a.attnum
                WHERE a.attrelid = %s::regclass AND a.attname = 'id'
            """, [partition_table])
            row = cursor.fetchone()
            self.assertIsNotNone(row, f"'{partition_table}.id' has no default")
            cursor.execute(f"SELECT {row[0]}")
            self.assertIsNotNone(cursor.fetchone()[0])

    def test_partition_id_default(self):
        self._assert_partition_id_auto_increments(self._make_gcc())

    def test_partition_id_default_with_non_standard_sequence_name(self):
        """ Postgres only calls it '<table>_id_seq' if that name was free when the column was created """
        with connection.cursor() as cursor:
            cursor.execute(f"ALTER SEQUENCE {self.BASE_TABLE_NAME}_id_seq RENAME TO {self.BASE_TABLE_NAME}_id_seq1")
        self._assert_partition_id_auto_increments(self._make_gcc())
