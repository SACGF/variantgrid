import logging
import os
import subprocess
import time
import uuid
from concurrent.futures.thread import ThreadPoolExecutor

import cyvcf2
import redis
from django.conf import settings
from django.core.cache import cache
from django.db import connection

from library.utils.file_utils import mk_path, name_from_filename
from snpdb.models import VariantCoordinate
from snpdb.variant_pk_lookup import VariantPKLookup
from snpdb.vcf_utils import vcf_get_ref_alt_svlen_and_modification
from upload.enums import UploadStepTaskType, VCFPipelineStage
from upload.models import ModifiedImportedVariant, UploadStep
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
                                            pipeline_stage=VCFPipelineStage.PRE_DATA_INSERTION,
                                            input_filename=unknown_variants_filename,
                                            items_to_process=len(unknown_variants))
    upload_step.launch_task(InsertUnknownVariantsTask)


def _handle_unknown_variants_in_thread(*args):
    """ The thread opens its own DB connection, which nothing else will close """
    try:
        handle_unknown_variants(*args)
    finally:
        connection.close()


class BulkUnknownVariantInserter:

    def __init__(self, upload_step):
        self.upload_step = upload_step
        self.unknown_dir = upload_step.upload_pipeline.get_pipeline_processing_subdir("vcf_unknown")
        self.rows_processed = 0
        self.unknown_variants_batch_id = 0
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.child_futures = []
        self.variant_pk_lookup = VariantPKLookup(upload_step.genome_build)

    def process_vcf_record(self, variant):
        # Pre-processed by vcf_filter_unknown_contigs so only recognised contigs present
        # This has been decomposed (only be 1 alt per line)
        ref, alt, svlen, _ = vcf_get_ref_alt_svlen_and_modification(variant,
                                                                    old_variant_info=ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
        variant_coordinate = VariantCoordinate(chrom=variant.CHROM, position=variant.POS, ref=ref, alt=alt, svlen=svlen)
        self.variant_pk_lookup.add(variant_coordinate)
        self.rows_processed += 1
        self.batch_process_check()

    def finish(self):
        self.batch_process_check(insert_all=True)

        # Make sure all child threads returned ok
        for f in self.child_futures:
            f.result()
        self.executor.shutdown()

    def batch_process_check(self, insert_all=False):
        if insert_all:
            minimum_size = 0
        else:
            minimum_size = settings.SQL_BATCH_INSERT_SIZE

        self.variant_pk_lookup.batch_check(minimum_size)

        # Unknown variants
        num_unknown_variants = len(self.variant_pk_lookup.unknown_variant_coordinates)
        if num_unknown_variants and num_unknown_variants >= minimum_size:
            uv_basename = f"unknown_variants_step_{self.upload_step.pk}_batch_{self.unknown_variants_batch_id}.csv"
            unknown_variants_filename = os.path.join(self.unknown_dir, uv_basename)
            huv_args = (self.upload_step.upload_pipeline, self.variant_pk_lookup.unknown_variant_coordinates,
                        unknown_variants_filename)

            if USE_THREADS:
                future = self.executor.submit(_handle_unknown_variants_in_thread, *huv_args)
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

        vcf_reader = cyvcf2.VCF(upload_step.input_filename)
        for i, v in enumerate(vcf_reader):
            try:
                bulk_inserter.process_vcf_record(v)
            except:
                logging.warning("SeparateUnknownVariantsTask died processing '%s' line %d record: %s",
                                upload_step.input_filename, i + 1, v)
                raise

        bulk_inserter.finish()  # Any leftovers
        return bulk_inserter.rows_processed


class AnnotateImportedVCFTask(ImportVCFStepTask):
    """ Indexes then annotates the VCF """
    def process_items(self, upload_step):
        cmd_index = ['bcftools', 'index', upload_step.input_filename]
        subprocess.run(cmd_index, capture_output=True, check=True)

        common_settings = settings.VCF_IMPORT_COMMON_FILTERS.get(upload_step.genome_build.name)
        gnomad_af_filename = common_settings["gnomad_af_filename"]
        if not gnomad_af_filename.startswith("/"):
            gnomad_af_filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, gnomad_af_filename)

        # Some VCFs come in with AF already populated, so we need to pick a unique INFO field
        column = f"INFO/{settings.VCF_IMPORT_COMMON_FILTER_INFO}:=INFO/AF"
        cmd_annotate = [
            "bcftools", "annotate", "-force",
            "--annotations", gnomad_af_filename, "--columns", column,
            "--output-type", "z", "--output", upload_step.output_filename,
            upload_step.input_filename
        ]
        subprocess.run(cmd_annotate, capture_output=True, check=True)
        return 1


class InsertUnknownVariantsTask(ImportVCFStepTask):
    """ Loci/Variants may have been inserted since the unknown variants were written, so we need to check
        in this process (which runs under a lock and in variant_id_single_workers so will avoid race conditions """
    LOCK_ID = "insert-unknown-variants-lock"
    LOCK_EXPIRE = 60 * 5  # 5 minutes
    # Longer than LOCK_EXPIRE, so a lock left by a killed worker expires and is taken, and long enough
    # for Redis to come back after a restart (BusyLoadingError while it reloads its dataset) #2189
    LOCK_WAIT = 60 * 15
    LOCK_POLL_SECONDS = 10

    def process_items(self, upload_step):
        if not settings.UPLOAD_ENABLED:
            raise ValueError("Uploads disabled, this should not have been called!")

        lock_token = self._acquire_lock()
        try:
            return InsertUnknownVariantsTask.process_items_in_lock(upload_step)
        finally:
            self._release_lock(lock_token)

    def _acquire_lock(self) -> str:
        """ Never skip the insert: the step would report success with its unknown variants missing.
            Returns the token to release with """
        lock_token = uuid.uuid4().hex
        deadline = time.monotonic() + self.LOCK_WAIT
        while True:
            try:
                if cache.add(self.LOCK_ID, lock_token, self.LOCK_EXPIRE):  # fails if already exists
                    return lock_token
                logging.info("%s is held, waiting", self.LOCK_ID)
            except redis.exceptions.ConnectionError as e:  # includes BusyLoadingError
                logging.warning("Redis unavailable acquiring %s: %s", self.LOCK_ID, e)
            if time.monotonic() > deadline:
                raise TimeoutError(f"Could not acquire '{self.LOCK_ID}' within {self.LOCK_WAIT} seconds")
            time.sleep(self.LOCK_POLL_SECONDS)

    def _release_lock(self, lock_token: str):
        """ Past LOCK_EXPIRE the lock may have been taken by another task - leave theirs alone """
        if cache.get(self.LOCK_ID) == lock_token:
            logging.info("Releasing %s", self.LOCK_ID)
            cache.delete(self.LOCK_ID)
        else:
            logging.warning("%s expired while held", self.LOCK_ID)

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
                chrom, position, ref, alt, svlen_str = line.strip().split(",")  # Not quoted, exactly 5 columns
                if svlen_str:
                    svlen = int(svlen_str)
                else:
                    svlen = None
                variant_coordinate = VariantCoordinate(chrom=chrom, position=int(position), ref=ref, alt=alt, svlen=svlen)
                variant_pk_lookup.add(variant_coordinate)
                variant_pk_lookup.batch_check(settings.SQL_BATCH_INSERT_SIZE, insert_unknown=True)
                items_processed += 1
            # Any remaining
            variant_pk_lookup.batch_check(insert_unknown=True)
        return items_processed


SeparateUnknownVariantsTask = app.register_task(SeparateUnknownVariantsTask())
AnnotateImportedVCFTask = app.register_task(AnnotateImportedVCFTask())
InsertUnknownVariantsTask = app.register_task(InsertUnknownVariantsTask())
