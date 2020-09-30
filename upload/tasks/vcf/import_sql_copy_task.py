from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.vcf import sql_copy_files
from variantgrid.celery import app


class ImportCohortGenotypeSQLCopyTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        table_name = upload_step.import_variant_table
        input_filename = upload_step.input_filename
        return sql_copy_files.cohort_genotype_sql_copy_csv(input_filename, table_name)


class ImportModifiedImportedVariantSQLCopyTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        input_filename = upload_step.input_filename
        return sql_copy_files.modified_imported_variant_sql_copy_csv(input_filename)


ImportCohortGenotypeSQLCopyTask = app.register_task(ImportCohortGenotypeSQLCopyTask())  # @UndefinedVariable
ImportModifiedImportedVariantSQLCopyTask = app.register_task(ImportModifiedImportedVariantSQLCopyTask())  # @UndefinedVariable
