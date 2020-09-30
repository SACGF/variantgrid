from django.conf import settings
from django.core.cache import cache
import logging
import os

from library.file_utils import name_from_filename, mk_path
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app


class InsertUnknownVariantsTask(ImportVCFStepTask):
    """ Loci/Variants may have been inserted since the unknown variants were written, so we need to check
        in this process (which runs under a lock and in variant_id_single_workers so will avoid race conditions """

    def process_items(self, upload_step):
        LOCK_EXPIRE = 60 * 5  # 5 minutes
        LOCK_ID = "insert-unknown-variants-lock"

        if not settings.UPLOAD_ENABLED:
            raise ValueError(f"Uploads disabled, this should not have been called!")

        if cache.add(LOCK_ID, "true", LOCK_EXPIRE):  # fails if already exists
            try:
                return InsertUnknownVariantsTask.process_items_in_lock(upload_step)
            finally:
                logging.info("Releasing %s", LOCK_ID)
                cache.delete(LOCK_ID)
        else:
            logging.info("Someone else is doing the import")

    @staticmethod
    def process_items_in_lock(upload_step):
        logging.info("InsertUnknownVariantsTask.process_items() - %s", upload_step)

        # TODO: Fix this, bit of a hack at the moment....
        input_filename = upload_step.input_filename
        logging.info("file: %s", input_filename)

        pipeline_processing_dir = upload_step.upload_pipeline.get_pipeline_processing_dir()
        working_dir = os.path.join(pipeline_processing_dir, "unknown_variants", name_from_filename(input_filename))
        mk_path(working_dir)

        variant_pk_lookup = VariantPKLookup(upload_step.genome_build, working_dir=working_dir)

        with open(input_filename) as f:
            items_processed = 0
            # Python CSV reader dies with extremely long lines, so we just do by hand (not quoted or anything)
            for line in f:
                chrom, position, ref, alt = line.strip().split(",")  # Not quoted, exactly 4 columns
                variant_pk_lookup.add(chrom, position, ref, alt)
                variant_pk_lookup.batch_check(settings.SQL_BATCH_INSERT_SIZE, insert_unknown=True)
                items_processed += 1
            # Any remaining
            variant_pk_lookup.batch_check(insert_unknown=True)
        return items_processed


InsertUnknownVariantsTask = app.register_task(InsertUnknownVariantsTask())  # @UndefinedVariable
