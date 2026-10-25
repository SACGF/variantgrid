import os
import tempfile

from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase

from library.utils.file_utils import IteratorFile
from snpdb.enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from upload.enums import ModifiedImportedVariantOperation, UploadedFileTypes
from upload.models.models import (
    FileUpload,
    ModifiedImportedVariant,
    SimpleVCFImportInfo,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.vcf.import_sql_copy_task import ImportCohortGenotypeSQLCopyTask
from upload.vcf import sql_copy_files
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


class ImportCohortGenotypeSQLCopyTaskTest(TestCase):
    """ A batch on the other side of a split boundary already copied variants 1 and 2 - the later rows are
        removed, and a later row with a different genotype is warned about on the upload """

    def setUp(self):
        user = User.objects.get_or_create(username="cg_copy_dedupe_user")[0]
        file_upload = FileUpload.objects.create(user=user, name="test.vcf", path="/tmp/test.vcf",
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.COMMAND_LINE)
        pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        UploadStep.objects.create(upload_pipeline=pipeline, name=UploadStep.PREPROCESS_VCF_NAME, sort_order=0)
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant_ids = [slowly_create_test_variant("3", 128198980 + i, 'A', 'T', grch37).pk for i in range(2)]
        with connection.cursor() as cursor:
            cursor.execute("CREATE TEMP TABLE cg_dedup_test (LIKE snpdb_cohortgenotype INCLUDING IDENTITY) "
                           "ON COMMIT DROP")
        earlier = self._write_csv([self._row(variant_id, "ER") for variant_id in self.variant_ids])
        sql_copy_files.cohort_genotype_sql_copy_csv(earlier, "cg_dedup_test")
        later = self._write_csv([self._row(self.variant_ids[0], "ER"), self._row(self.variant_ids[1], "OE")])
        self.upload_step = UploadStep.objects.create(upload_pipeline=pipeline, name="copy", sort_order=1,
                                                     import_variant_table="cg_dedup_test", input_filename=later)

    @staticmethod
    def _row(variant_id, samples_zygosity: str) -> str:
        # COPY text format: arrays escape their commas
        ref, het, hom = (samples_zygosity.count(z) for z in "REO")
        array = r"{20\,30}"
        return f"7,{variant_id},,{ref},{het},{hom},0,{samples_zygosity},{array},{array},{array},{array},{array},,[],{{}}\n"

    def _write_csv(self, rows) -> str:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.writelines(rows)
            self.addCleanup(os.unlink, f.name)
            return f.name

    def test_conflicting_duplicate_warns(self):
        ImportCohortGenotypeSQLCopyTask.process_items(self.upload_step)

        with connection.cursor() as cursor:
            cursor.execute("SELECT variant_id, samples_zygosity FROM cg_dedup_test ORDER BY id")
            self.assertEqual(cursor.fetchall(), [(self.variant_ids[0], "ER"), (self.variant_ids[1], "ER")])
        rmdup = ModifiedImportedVariant.objects.filter(operation=ModifiedImportedVariantOperation.RMDUP)
        self.assertEqual({(self.variant_ids[0], None),
                          (self.variant_ids[1], ImportCohortGenotypeSQLCopyTask.CONFLICTING_GENOTYPES_DETAIL)},
                         set(rmdup.values_list("variant_id", "operation_detail")))
        info = SimpleVCFImportInfo.objects.get(type=SimpleVCFImportInfo.CONFLICTING_GENOTYPES)
        self.assertEqual(1, info.count)
