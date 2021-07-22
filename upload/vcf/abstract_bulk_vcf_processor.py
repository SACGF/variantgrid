from django.conf import settings
import os

from library.django_utils.django_file_utils import get_import_processing_filename
from snpdb.models import Variant
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadStep, ModifiedImportedVariant, UploadStepTaskType, VCFPipelineStage
from upload.tasks.vcf.import_sql_copy_task import ImportModifiedImportedVariantSQLCopyTask
from upload.vcf.sql_copy_files import write_sql_copy_csv


class AbstractBulkVCFProcessor:
    """ The minimum a VCF processor needs to do is
        set max_variant_id and insert multi-allelics so we know what was normalised """

    def __init__(self, upload_step, preprocess_vcf_import_info, batch_size=settings.SQL_BATCH_INSERT_SIZE):
        self.upload_step = upload_step
        self.upload_pipeline = upload_step.upload_pipeline
        self.preprocess_vcf_import_info = preprocess_vcf_import_info
        self.batch_size = batch_size

        self.rows_processed = 0
        self.max_variant_id = 0
        self.modified_imported_variants_file_id = 0
        self.variant_hashes = []
        self.modified_imported_variant_hashes = []
        self.modified_imported_variants = []
        self.variant_pk_lookup = VariantPKLookup.factory(upload_step.genome_build)

    @property
    def genome_build(self):
        return self.upload_pipeline.genome_build

    def set_max_variant(self, variant_ids):
        max_returned_variant_id = max(map(int, variant_ids))
        self.max_variant_id = max(self.max_variant_id, max_returned_variant_id)

    def get_max_variant_id(self):
        """ 0 means it was never set, so we return None """
        return self.max_variant_id or None

    def add_modified_imported_variant(self, variant, variant_hash, miv_hash_list=None, miv_list=None):
        old_multiallelic = variant.INFO.get("OLD_MULTIALLELIC")
        old_variant = variant.INFO.get("OLD_VARIANT")

        if old_multiallelic or old_variant:
            if miv_hash_list is None:
                miv_hash_list = self.modified_imported_variant_hashes
            if miv_list is None:
                miv_list = self.modified_imported_variants

            miv_hash_list.append(variant_hash)
            if old_variant:
                old_variant_formatted = ModifiedImportedVariant.format_old_variant(old_variant, self.genome_build)
            else:
                old_variant_formatted = [None]
            for ov in old_variant_formatted:
                miv_list.append((old_multiallelic, old_variant, ov))

    def process_modified_imported_variants(self, variant_ids_by_hash):
        modified_imported_variants = []
        for variant_hash, modified_imported_variants_tuple in zip(self.modified_imported_variant_hashes, self.modified_imported_variants):
            variant_id = variant_ids_by_hash[variant_hash]
            modified_imported_variants.append((self.preprocess_vcf_import_info.pk, variant_id) + modified_imported_variants_tuple)

        modified_imported_variants_basename = f"modified_imported_variants_step_{self.upload_step.pk}_batch_{self.modified_imported_variants_file_id}.csv"
        modified_imported_variants_filename = get_import_processing_filename(self.upload_pipeline.pk, modified_imported_variants_basename)
        write_sql_copy_csv(modified_imported_variants, modified_imported_variants_filename)
        self.create_modified_imported_variants_job(len(modified_imported_variants), modified_imported_variants_filename)

        self.modified_imported_variants_file_id += 1
        self.modified_imported_variant_hashes = []
        self.modified_imported_variants = []

    def create_modified_imported_variants_job(self, num_modified_imported_variants, input_filename):
        # Should this be it's own table?? Should only be 10s of thousands per file - and we don't need it fast?

        if not os.path.exists(input_filename):
            msg = f"create_modified_imported_variants_job: input file: '{input_filename}' does not exist."
            raise ValueError(msg)

        name = "ModifiedImportedVariant SQL COPY"
        sort_order = self.upload_pipeline.get_max_step_sort_order() + 1
        sql_job = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                            name=name,
                                            sort_order=sort_order,
                                            task_type=UploadStepTaskType.SQL,
                                            pipeline_stage=VCFPipelineStage.DATA_INSERTION,
                                            input_filename=input_filename,
                                            items_to_process=num_modified_imported_variants)

        sql_job.launch_task(ImportModifiedImportedVariantSQLCopyTask)

    @staticmethod
    def get_ref_alt(variant):
        ref = variant.REF.strip().upper()
        if variant.ALT:
            alt = variant.ALT[0].strip().upper()
        else:
            alt = Variant.REFERENCE_ALT
        return ref, alt
