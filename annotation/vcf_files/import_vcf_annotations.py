#!/usr/bin/env python3
import logging
import os

from django.conf import settings
import cyvcf2
from django.utils import timezone

from annotation.annotation_version_querysets import get_unannotated_variants_qs
from annotation.models import AnnotationRun, re, VEPSkippedReason
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from library.vcf_utils import cyvcf2_header_types


def import_vcf_annotations(annotation_run: AnnotationRun, insert_variants=True,
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
        annotation_run.annotated_count = bulk_inserter.rows_processed

        handle_vep_warnings(annotation_run, bulk_inserter)

        if delete_temp_files:
            bulk_inserter.remove_processing_files()

        annotation_run.upload_end = timezone.now()
        annotation_run.save()
        logging.info("%d variants imported!", annotation_run.annotated_count)
    except:
        annotation_run.set_task_log("BulkVEPVCFAnnotationInserter line number", bulk_inserter.rows_processed)
        raise


def handle_vep_warnings(annotation_run: AnnotationRun, bulk_inserter):
    """ set annotation_run fields 'vep_warnings' and 'vep_skipped_count'
        Create VariantAnnotation entries for any skipped variants, filling in vep_skipped_reason

        Uses bulk_inserter as that handles DB partitions """

    annotation_run.vep_skipped_count = annotation_run.annotation_range_lock.count - annotation_run.annotated_count
    if annotation_run.vep_skipped_count:
        vep_warning_filename = annotation_run.vcf_annotated_filename + "_warnings.txt"
        print(f"Looking for '{vep_warning_filename}'")
        if os.path.exists(vep_warning_filename):
            with open(vep_warning_filename) as f:
                annotation_run.vep_warnings = f.read()

        # Process warnings file, and work out reasons (if any)
        incomplete_variant_ids = set()
        skipped_contigs = set()
        if annotation_run.vep_warnings:
            next_line_incomplete = False
            for line in annotation_run.vep_warnings.splitlines():
                if m := re.match(r"WARNING: Chromosome (\S+) not found", line):
                    skipped_contigs.add(m.group(1))
                elif re.match(r"WARNING: VCF line on line \d+ looks incomplete, skipping:", line):
                    next_line_incomplete = True
                elif next_line_incomplete:
                    if m := re.match(r".*variant_id=(\d+)"):
                        incomplete_variant_ids.add(int(m.group(1)))
                    next_line_incomplete = False

        print(f"skipped_contigs = {skipped_contigs}")

        version = annotation_run.annotation_range_lock.version
        annotation_version = version.get_any_annotation_version()
        for v in get_unannotated_variants_qs(annotation_version,
                                             min_variant_id=annotation_run.annotation_range_lock.min_variant_id,
                                             max_variant_id=annotation_run.annotation_range_lock.max_variant_id):
            if v.locus.contig.name in skipped_contigs:
                reason = VEPSkippedReason.UNKNOWN_CONTIG
            elif v.pk in incomplete_variant_ids:
                reason = VEPSkippedReason.INCOMPLETE
            else:
                reason = VEPSkippedReason.UNKNOWN

            va = {"version_id": version.pk, "variant_id": v.pk,
                  "annotation_run_id": annotation_run.pk, "vep_skipped_reason": reason,
                  "predictions_num_pathogenic": 0, "predictions_num_benign": 0}
            bulk_inserter.variant_annotation_list.append(va)

        # Only insert the columns we care about
        bulk_inserter.all_variant_columns = ["vep_skipped_reason"]
        bulk_inserter.finish()
