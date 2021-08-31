import gzip
import logging
import os
from concurrent.futures.thread import ThreadPoolExecutor

from django.conf import settings
from django.core.cache import cache

from library.file_utils import name_from_filename, mk_path
from snpdb import variant_collection
from snpdb.models import Variant
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadStep, UploadStepTaskType, VCFPipelineStage, VCFSkippedGVCFNonVarBlocks
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.vcf import sql_copy_files
from variantgrid.celery import app

USE_THREADS = True


def handle_unknown_variants(upload_pipeline, unknown_variants, unknown_variants_filename):
    sql_copy_files.write_sql_copy_csv(unknown_variants, unknown_variants_filename)

    name = UploadStep.CREATE_UNKNOWN_LOCI_AND_VARIANTS_TASK_NAME
    sort_order = upload_pipeline.get_max_step_sort_order() + 1
    upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                            name=name,
                                            sort_order=sort_order,
                                            task_type=UploadStepTaskType.CELERY,
                                            pipeline_stage=VCFPipelineStage.INSERT_UNKNOWN_VARIANTS,
                                            input_filename=unknown_variants_filename,
                                            items_to_process=len(unknown_variants))
    upload_step.launch_task(InsertUnknownVariantsTask)


class BulkUnknownVariantInserter:

    def __init__(self, upload_step):
        self.upload_step = upload_step
        self.unknown_dir = upload_step.upload_pipeline.get_pipeline_processing_subdir("vcf_unknown")
        self.rows_processed = 0
        self.unknown_variants_batch_id = 0
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.child_futures = []
        self.variant_pk_lookup = VariantPKLookup(upload_step.genome_build)
        # Settings
        self.store_gvcf_non_var_blocks = settings.VCF_IMPORT_STORE_GVCF_NON_VAR_BLOCKS
        self.num_skipped_gvcf_non_var_blocks = 0

    def process_vcf_line(self, line):
        columns = line.split("\t")
        if len(columns) < 5:
            msg = f"VCF line: '{line}' has < 5 columns"
            raise ValueError(msg)

        # Pre-processed by vcf_filter_unknown_contigs so only recognised contigs present
        chrom = columns[0]
        position = columns[1]
        ref = columns[3]
        alt = columns[4]
        if self.store_gvcf_non_var_blocks is False:
            if alt == "<NON_REF>":
                self.num_skipped_gvcf_non_var_blocks += 1
                return  # Skip

        if Variant.is_ref_alt_reference(ref, alt):
            alt = Variant.REFERENCE_ALT

        # Always insert a reference with every non-ref alt variant
        for alt in {alt, Variant.REFERENCE_ALT}:
            self.variant_pk_lookup.add(chrom, position, ref, alt)
        self.batch_process_check()

    def finish(self):
        self.batch_process_check(insert_all=True)

        # Make sure all child threads returned ok
        for f in self.child_futures:
            f.result()

        if self.num_skipped_gvcf_non_var_blocks:
            VCFSkippedGVCFNonVarBlocks.objects.create(upload_step=self.upload_step,
                                                      num_skipped=self.num_skipped_gvcf_non_var_blocks)

    def batch_process_check(self, insert_all=False):
        if insert_all:
            minimum_redis_pipeline_size = 0
            minimum_unknown_insert_size = 0
        else:
            minimum_redis_pipeline_size = settings.REDIS_PIPELINE_SIZE
            minimum_unknown_insert_size = settings.SQL_BATCH_INSERT_SIZE

        self.variant_pk_lookup.batch_check(minimum_redis_pipeline_size)

        # Unknown variants
        num_unknown_variants = len(self.variant_pk_lookup.unknown_variant_coordinates)
        if num_unknown_variants and num_unknown_variants >= minimum_unknown_insert_size:
            unknown_variants_basename = f"unknown_variants_step_{self.upload_step.pk}_batch_{self.unknown_variants_batch_id}.csv"
            unknown_variants_filename = os.path.join(self.unknown_dir, unknown_variants_basename)
            huv_args = (self.upload_step.upload_pipeline, self.variant_pk_lookup.unknown_variant_coordinates, unknown_variants_filename)

            if USE_THREADS:
                future = self.executor.submit(handle_unknown_variants, *huv_args)
                self.child_futures.append(future)
            else:
                handle_unknown_variants(*huv_args)

            self.variant_pk_lookup.clear()
            self.unknown_variants_batch_id += 1


class SeparateUnknownVariantsTask(ImportVCFStepTask):
    """ This is run on split VCF file.
        If it finds unknown, it writes files and launches InsertUnknownVariantsTask on them """

    def process_items(self, upload_step):
        bulk_inserter = BulkUnknownVariantInserter(upload_step)
        with gzip.open(upload_step.input_filename, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                bulk_inserter.process_vcf_line(line)

            bulk_inserter.finish()  # Any leftovers

        variant_collection.count = bulk_inserter.rows_processed
        return bulk_inserter.rows_processed


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


SeparateUnknownVariantsTask = app.register_task(SeparateUnknownVariantsTask())
InsertUnknownVariantsTask = app.register_task(InsertUnknownVariantsTask())
