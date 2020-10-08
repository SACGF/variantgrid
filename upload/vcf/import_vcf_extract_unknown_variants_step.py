import gzip
from concurrent.futures.thread import ThreadPoolExecutor
from django.conf import settings
import itertools
import logging
import os

from snpdb import variant_collection
from snpdb.models import Variant
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadStep, UploadStepTaskType, \
    VCFPipelineStage, UploadStepMultiFileOutput, VCFSkippedGVCFNonVarBlocks
from upload.tasks.vcf.unknown_variants_task import InsertUnknownVariantsTask
from upload.vcf import sql_copy_files

USE_THREADS = True


def write_split_vcf(upload_step, vcf_header_lines, vcf_lines, split_vcf_filename, items_to_process):
    with gzip.open(split_vcf_filename, "wt") as f:
        f.writelines(vcf_header_lines)
        f.writelines(vcf_lines)

    UploadStepMultiFileOutput.objects.create(upload_step=upload_step,
                                             output_filename=split_vcf_filename,
                                             items_to_process=items_to_process)


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

    def __init__(self, upload_step, vcf_header_lines):
        self.upload_step = upload_step
        self.unknown_dir = upload_step.upload_pipeline.get_pipeline_processing_subdir("vcf_unknown")
        self.split_vcf_dir = upload_step.upload_pipeline.get_pipeline_processing_subdir("split_vcf")
        self.rows_processed = 0
        self.vcf_header_lines = vcf_header_lines
        self.vcf_lines = []
        self.vcf_locus_lines = []  # Locus buffers used to avoid splitting locus across files
        self.last_locus_tuple = None
        self.split_vcf_batch_id = 0
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

        # Need to avoid splitting locus across multiple files
        locus_tuple = (chrom, position, ref)
        self.vcf_locus_lines.append(line)
        if self.last_locus_tuple:
            if self.last_locus_tuple != locus_tuple:
                self.finish_locus()
        self.last_locus_tuple = locus_tuple

        if alt == ref:
            alt = Variant.REFERENCE_ALT

        # Always insert a reference with every non-ref alt variant
        for alt in {alt, Variant.REFERENCE_ALT}:
            self.variant_pk_lookup.add(chrom, position, ref, alt)
        self.batch_process_check()

    def finish_locus(self):
        self.vcf_lines.extend(self.vcf_locus_lines)
        self.vcf_locus_lines = []

    def finish(self):
        self.finish_locus()
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
            minimum_split_file_size = 0
        else:
            minimum_redis_pipeline_size = settings.REDIS_PIPELINE_SIZE
            minimum_unknown_insert_size = settings.SQL_BATCH_INSERT_SIZE
            minimum_split_file_size = settings.VCF_IMPORT_FILE_SPLIT_ROWS

        # These operations are I/O heavy (writing files, calling redis) so use threads
        num_vcf_lines = len(self.vcf_lines)
        if num_vcf_lines and num_vcf_lines >= minimum_split_file_size:
            split_vcf_basename = f"split_vcf_step_{self.upload_step.pk}_batch_{self.split_vcf_batch_id}.vcf.gz"
            split_vcf_filename = os.path.join(self.split_vcf_dir, split_vcf_basename)
            wpv_args = (self.upload_step, self.vcf_header_lines, self.vcf_lines, split_vcf_filename, num_vcf_lines)
            if USE_THREADS:
                future = self.executor.submit(write_split_vcf, *wpv_args)
                self.child_futures.append(future)
            else:
                write_split_vcf(*wpv_args)

            self.split_vcf_batch_id += 1
            self.vcf_lines = []

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


def import_vcf_extract_unknown_variants_and_split_file(upload_step, vcf_file):
    """ As quickly as possible - go through VCF and find any unknown variants (then we can insert them)
        Splits VCF up into multiple output files, so they can be processed in parallel
        returns items_processed
    """

    logging.debug("vcf_separate_known_unknown_variants start")
    vcf_header_lines = []
    line = None
    for line in vcf_file:
        if line.startswith("#"):
            vcf_header_lines.append(line)
            line = None  # So it won't be processed if nothing but header
        else:
            break

    bulk_inserter = BulkUnknownVariantInserter(upload_step, vcf_header_lines)
    if line:
        for line in itertools.chain([line], vcf_file):
            bulk_inserter.process_vcf_line(line)

        bulk_inserter.finish()  # Any leftovers

    variant_collection.count = bulk_inserter.rows_processed
    return bulk_inserter.rows_processed
