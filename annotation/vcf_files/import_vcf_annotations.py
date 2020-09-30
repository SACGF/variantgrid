#!/usr/bin/env python3
import logging

from django.conf import settings
import cyvcf2
from django.utils import timezone

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from library.log_utils import get_traceback
from library.vcf_utils import cyvcf2_header_types


def import_vcf_annotations(annotation_run, insert_variants=True,
                           delete_temp_files=settings.VARIANT_ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS):
    annotation_run.upload_start = timezone.now()
    annotation_run.upload_attempts += 1
    annotation_run.save()

    vcf_reader = cyvcf2.VCF(annotation_run.vcf_annotated_filename)
    header_types = cyvcf2_header_types(vcf_reader)
    infos = header_types["INFO"]
    bulk_inserter = BulkVEPVCFAnnotationInserter(annotation_run, infos, insert_variants)
    if annotation_run.upload_attempts > 1:
        logging.warning(f"Upload attempt {annotation_run.upload_attempts} - deleting data")
        bulk_inserter.remove_processing_files()

    try:
        for v in vcf_reader:
            bulk_inserter.process_entry(v)

        bulk_inserter.finish()  # Any leftovers

        if delete_temp_files:
            bulk_inserter.remove_processing_files()

        annotation_run.variant_count = bulk_inserter.rows_processed
        annotation_run.upload_end = timezone.now()
        annotation_run.save()
        logging.info("%d variants imported!", annotation_run.variant_count)
    except:
        annotation_run.set_task_log("BulkVEPVCFAnnotationInserter line number", bulk_inserter.rows_processed)
        raise
