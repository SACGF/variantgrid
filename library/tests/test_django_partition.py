"""
Tests for RelatedModelsPartitionModel partition creation and deletion.
"""

from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from genes.models import GeneCoverage, GeneCoverageCollection
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


class PartitionDeleteTests(TestCase):
    """ Records are deleted by dropping the partition, so the related models point at the collection with
        DO_NOTHING. A CASCADE there makes Django's collector issue a DELETE against the inheritance parent,
        which locks every other collection's partition and deadlocks concurrent imports
        (SACGF/variantgrid_sapath#450). """

    BASE_TABLE_NAME = "genes_genecoverage"

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch38()

    def _make_gcc(self) -> GeneCoverageCollection:
        return GeneCoverageCollection.objects.create(path="/tmp/nonexistent.tsv",
                                                     data_state=DataState.COMPLETE,
                                                     genome_build=self.genome_build)

    def _make_base_table_record(self, gcc: GeneCoverageCollection) -> GeneCoverage:
        """ The model's db_table is the inheritance parent, so a plain create bypasses the partition """
        return GeneCoverage.objects.create(gene_coverage_collection=gcc, original_gene_symbol="RUNX1",
                                           original_transcript="NM_001754.4", min=0, mean=0.0, std_dev=0.0)

    def _count_base_table_records(self, gcc: GeneCoverageCollection) -> int:
        with connection.cursor() as cursor:
            cursor.execute(f"SELECT count(*) FROM ONLY {self.BASE_TABLE_NAME} "
                           "WHERE gene_coverage_collection_id = %s", [gcc.pk])
            return cursor.fetchone()[0]

    def test_delete_does_not_delete_across_inheritance_parent(self):
        gcc = self._make_gcc()
        with CaptureQueriesContext(connection) as queries:
            gcc.delete()

        base_table = f'"{self.BASE_TABLE_NAME}"'  # quoted, so it doesn't also match the collection table
        tree_wide = [q["sql"] for q in queries
                     if q["sql"].lstrip().upper().startswith("DELETE")
                     and base_table in q["sql"] and f"ONLY {base_table}" not in q["sql"]]
        self.assertEqual(tree_wide, [], f"Delete expanded over the inheritance tree: {tree_wide}")

    def test_delete_removes_records_left_in_base_table(self):
        gcc = self._make_gcc()
        self._make_base_table_record(gcc)
        self.assertEqual(self._count_base_table_records(gcc), 1)

        gcc.delete()
        self.assertEqual(self._count_base_table_records(gcc), 0)
