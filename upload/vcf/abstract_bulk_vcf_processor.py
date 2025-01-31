import abc
import os
from collections import Counter
from typing import Optional

import cyvcf2
from django.conf import settings

from library.django_utils.django_file_utils import get_import_processing_filename
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.genomics.vcf_utils import vcf_get_ref_alt_svlen_and_modification
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadStep, ModifiedImportedVariant, UploadStepTaskType, VCFPipelineStage, \
    SimpleVCFImportInfo, ModifiedImportedVariantOperation
from upload.tasks.vcf.import_sql_copy_task import ImportModifiedImportedVariantSQLCopyTask
from upload.vcf.sql_copy_files import write_sql_copy_csv


class AbstractBulkVCFProcessor(abc.ABC):
    """ The minimum a VCF processor needs to do is
        set max_variant_id and insert multi-allelics so we know what was normalised """

    def __init__(self, upload_step, preprocess_vcf_import_info, batch_size=settings.SQL_BATCH_INSERT_SIZE):
        self.upload_step = upload_step
        self.upload_pipeline = upload_step.upload_pipeline
        self.preprocess_vcf_import_info = preprocess_vcf_import_info
        self.batch_size = batch_size

        self.rows_processed = 0
        self.max_variant_id = 0  # Non reference variant (ie that will be annotated)
        self.set_max_variant_called = False
        self.modified_imported_variants_file_id = 0
        self.variant_hashes = []
        self.modified_imported_variant_hashes = []
        self.modified_imported_variants = []
        self.variant_pk_lookup = VariantPKLookup(upload_step.genome_build)
        self.svlen_modifications = Counter()

    @abc.abstractmethod
    def process_entry(self, variant: cyvcf2.Variant):
        pass

    @abc.abstractmethod
    def _finish(self):
        pass

    def finish(self):
        self._finish()
        if self.svlen_modifications:
            for modification, count in self.svlen_modifications.items():
                message_string = f"{modification}: {count}"
                SimpleVCFImportInfo.objects.create(type=SimpleVCFImportInfo.SVLEN_MODIFIED, has_more_details=True,
                                                   upload_step=self.upload_step, message_string=message_string)

    @property
    def genome_build(self):
        return self.upload_pipeline.genome_build

    def get_ref_alt_svlen(self, variant) -> tuple[str, str, Optional[int]]:
        """ Ensures SVLEN fits symbolic alt, some callers eg Manta (non-Dragen) write positive SVLEN for dels
            We don't do this in vcf_clean_and_filter as we don't pull apart the INFO there
        """
        ref, alt, svlen, modification = vcf_get_ref_alt_svlen_and_modification(variant)
        if modification:
            self.svlen_modifications[modification] += 1
        return ref, alt, svlen

    def set_max_variant(self, variant_hashes, variant_ids):
        # Keep track of max annotated variant (only non-reference are annotated)
        self.set_max_variant_called = True
        non_ref_variant_ids = self.variant_pk_lookup.filter_non_reference(variant_hashes, variant_ids)
        if non_ref_variant_ids:
            max_returned_variant_id = max(map(int, non_ref_variant_ids))
            self.max_variant_id = max(self.max_variant_id, max_returned_variant_id)

    def get_max_variant_id(self):
        """ 0 means it was never set, so we return None """

        if self.max_variant_id:
            return self.max_variant_id

        if self.rows_processed and self.set_max_variant_called is False:
            msg = "max_variant_id not set, importers must call set_max_variant() if you process any variants!"
            raise ValueError(msg)

        return None  # No variants, or only reference

    def add_modified_imported_variant(self, variant: cyvcf2.Variant, variant_hash, miv_hash_list=None, miv_list=None):
        # This used to handle VT tags: OLD_MULTIALLELIC / OLD_VARIANT but now we handle BCFTOOLS only
        if bcftools_old_variant := variant.INFO.get(ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG):
            svlen = variant.INFO.get("SVLEN")

            if miv_hash_list is None:
                miv_hash_list = self.modified_imported_variant_hashes
            if miv_list is None:
                miv_list = self.modified_imported_variants

            if "," in bcftools_old_variant:  # Was a multi-allelic
                old_multiallelic = bcftools_old_variant
            else:
                old_multiallelic = None

            old_position = int(bcftools_old_variant.split("|")[1])
            if old_position != variant.POS:
                old_variant = bcftools_old_variant
            else:
                old_variant = None  # Wasn't normalized

            for ov in ModifiedImportedVariant.bcftools_format_old_variant(bcftools_old_variant, svlen, self.genome_build):
                # These 2 need to be in sync
                miv_hash_list.append(variant_hash)
                miv_list.append((ModifiedImportedVariantOperation.NORMALIZATION, old_multiallelic, old_variant, ov))

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
