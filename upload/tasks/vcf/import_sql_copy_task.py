from upload.enums import ModifiedImportedVariantOperation
from upload.models.models import (
    ModifiedImportedVariant,
    ModifiedImportedVariants,
    SimpleVCFImportInfo,
)
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.vcf import sql_copy_files
from variantgrid.celery import app


class ImportCohortGenotypeSQLCopyTask(ImportVCFStepTask):
    CONFLICTING_GENOTYPES_DETAIL = "conflicting genotype"
    CONFLICTING_GENOTYPES_MESSAGE = ("Same variant had conflicting genotypes in this VCF (an import / normalisation "
                                     "bug) - kept the first, removed the later")

    def process_items(self, upload_step):
        table_name = upload_step.import_variant_table
        input_filename = upload_step.input_filename
        num_copied = sql_copy_files.cohort_genotype_sql_copy_csv(input_filename, table_name)
        if removed := sql_copy_files.delete_duplicate_cohort_genotypes(input_filename, table_name):
            import_info = ModifiedImportedVariants.get_for_pipeline(upload_step.upload_pipeline)
            ModifiedImportedVariant.objects.bulk_create([
                ModifiedImportedVariant(import_info=import_info, variant_id=variant_id,
                                        operation=ModifiedImportedVariantOperation.RMDUP,
                                        operation_detail=self.CONFLICTING_GENOTYPES_DETAIL if conflicting else None)
                for variant_id, conflicting in removed
            ])
            if num_conflicting := sum(conflicting for _, conflicting in removed):
                SimpleVCFImportInfo.add_message_count(num_conflicting, self.CONFLICTING_GENOTYPES_MESSAGE, upload_step,
                                                      type=SimpleVCFImportInfo.CONFLICTING_GENOTYPES)
        return num_copied


class ImportModifiedImportedVariantSQLCopyTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        input_filename = upload_step.input_filename
        return sql_copy_files.modified_imported_variant_sql_copy_csv(input_filename)


ImportCohortGenotypeSQLCopyTask = app.register_task(ImportCohortGenotypeSQLCopyTask())  # @UndefinedVariable
ImportModifiedImportedVariantSQLCopyTask = app.register_task(ImportModifiedImportedVariantSQLCopyTask())  # @UndefinedVariable
