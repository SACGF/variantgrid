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

    def process_items(self, upload_step):
        table_name = upload_step.import_variant_table
        input_filename = upload_step.input_filename
        num_copied = sql_copy_files.cohort_genotype_sql_copy_csv(input_filename, table_name)
        if removed := sql_copy_files.delete_duplicate_cohort_genotypes(input_filename, table_name):
            import_info = ModifiedImportedVariants.get_for_pipeline(upload_step.upload_pipeline)
            modified_imported_variants = []
            num_conflicting = 0
            for variant_id, conflict_detail in removed:
                if conflict_detail:
                    operation = ModifiedImportedVariantOperation.CONFLICTING_DUPLICATE
                    num_conflicting += 1
                else:
                    operation = ModifiedImportedVariantOperation.RMDUP
                modified_imported_variants.append(
                    ModifiedImportedVariant(import_info=import_info, variant_id=variant_id, operation=operation,
                                            operation_detail=conflict_detail))
            ModifiedImportedVariant.objects.bulk_create(modified_imported_variants)
            if num_conflicting:
                SimpleVCFImportInfo.add_conflicting_genotypes(num_conflicting, upload_step)
        return num_copied


class ImportModifiedImportedVariantSQLCopyTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        input_filename = upload_step.input_filename
        return sql_copy_files.modified_imported_variant_sql_copy_csv(input_filename)


ImportCohortGenotypeSQLCopyTask = app.register_task(ImportCohortGenotypeSQLCopyTask())  # @UndefinedVariable
ImportModifiedImportedVariantSQLCopyTask = app.register_task(ImportModifiedImportedVariantSQLCopyTask())  # @UndefinedVariable
