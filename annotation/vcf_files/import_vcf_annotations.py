#!/usr/bin/env python3
import logging
import os
from typing import Optional

from django.conf import settings
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.models import AnnotationRun, VEPSkippedReason
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.vcf_files.bulk_annotsv_tsv_inserter import import_annotsv_tsv
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vep_annotation import vep_check_annotated_file_version_match
from annotation.vep_warnings import (
    TruncatedVEPOutputError,
    parse_skipped_variants_file,
    parse_vep_warnings,
)
from library.log_utils import log_traceback
from snpdb.models import VariantCoordinate


def import_vcf_annotations(
        annotation_run: AnnotationRun,
        insert_variants: bool = True,
        vep_version_check: bool = True,
        delete_temp_files: Optional[bool] = None):
    import cyvcf2

    from library.genomics.vcf_utils import cyvcf2_header_types

    if delete_temp_files is None:
        delete_temp_files = settings.IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS

    annotation_run.upload_start = timezone.now()
    annotation_run.upload_attempts += 1
    annotation_run.save()

    if vep_version_check:
        vep_check_annotated_file_version_match(annotation_run.variant_annotation_version, annotation_run.vcf_annotated_filename)
    vcf_reader = cyvcf2.VCF(annotation_run.vcf_annotated_filename)
    header_types = cyvcf2_header_types(vcf_reader)
    infos = header_types["INFO"]
    bulk_inserter = BulkVEPVCFAnnotationInserter(annotation_run, infos, insert_variants)
    if annotation_run.upload_attempts > 1:
        # Genuine retry of an existing run - clear any scratch files left by the previous attempt so
        # write_sql_copy_csv doesn't abort on "We don't want to overwrite ...". We deliberately only do
        # this on a real retry (upload_attempts > 1) - a first attempt that finds existing scratch files
        # means the import-processing dir is out of sync with the DB (eg moved dump / restored backup),
        # which the overwrite guard is meant to surface rather than silently delete. See #1596.
        logging.warning(f"Upload attempt {annotation_run.upload_attempts} - deleting data")
        bulk_inserter.remove_processing_files()

    try:
        # Counted separately from bulk_inserter.rows_processed, which only counts records carrying a CSQ
        records_read = 0
        try:
            for v in vcf_reader:
                records_read += 1
                bulk_inserter.process_entry(v)
        except Exception as e:
            # #1701: a cut mid-record surfaces from cyvcf2 as a bare htslib error code. Say what it is.
            msg = (f"AnnotationRun {annotation_run.pk}: failed reading "
                   f"'{annotation_run.vcf_annotated_filename}' after {records_read} records - "
                   f"the annotated VCF looks truncated. Retry the full run (not upload only). Cause: {e}")
            raise TruncatedVEPOutputError(msg) from e

        bulk_inserter.finish()  # Any leftovers
        annotation_run.annotated_count = bulk_inserter.rows_processed

        # Before handle_vep_skipped, which would otherwise write a skip row for every lost record
        check_all_dumped_records_accounted_for(annotation_run, records_read)

        handle_vep_skipped(annotation_run, bulk_inserter)

        if delete_temp_files:
            bulk_inserter.remove_processing_files()

        if annotation_run.annotsv_tsv_filename:
            try:
                import_annotsv_tsv(annotation_run)
            except Exception as e:
                # Best-effort: AnnotSV import errors must not fail the whole run.
                log_traceback()
                annotation_run.annotsv_error = (str(e) or "")[:100_000]
                annotation_run.annotsv_imported = False
                annotation_run.save()

        annotation_run.upload_end = timezone.now()
        annotation_run.save()
        logging.info("%d variants imported!", annotation_run.annotated_count)
    except:
        annotation_run.set_task_log("BulkVEPVCFAnnotationInserter line number", bulk_inserter.rows_processed)
        raise


def _vep_skipped_count(annotation_run: AnnotationRun) -> int:
    """ Records this run handed VEP that came back without annotation - its own arithmetic.

        annotation_range_lock.count is the fallback for legacy rows written before dump_count existed:
        it counts unannotated variants of *every* pipeline type in the range, so it over-reports for a
        run of one pipeline type over a mixed range. """
    annotated_count = annotation_run.annotated_count or 0
    if annotation_run.dump_count is not None:
        return annotation_run.dump_count - annotated_count
    return annotation_run.annotation_range_lock.count - annotated_count


def load_vep_warnings(annotation_run: AnnotationRun) -> Optional[str]:
    """ Read the `_warnings.txt` VEP writes next to the annotated VCF onto the run, returning it.

        Both the truncation check and handle_vep_skipped need it, and the check wants it on the row even
        when it fails - a warnings blob is the first thing you want on an ERROR run. Read once per import:
        the check runs first, so by the time handle_vep_skipped asks, the run already carries it. """
    if annotation_run.vep_warnings is None and annotation_run.vcf_annotated_filename:
        vep_warning_filename = os.path.abspath(annotation_run.vcf_annotated_filename + "_warnings.txt")
        logging.info("Looking for '%s'", vep_warning_filename)
        if os.path.exists(vep_warning_filename):
            with open(vep_warning_filename) as f:
                annotation_run.vep_warnings = f.read()
    return annotation_run.vep_warnings


def count_vep_skipped(annotation_run: AnnotationRun) -> tuple[int, str]:
    """ How many records VEP told us it dropped, and where that number came from.

        Prefers the structured --skipped_variants_file (#1701) - an exact row per skip, no pattern
        matching - and falls back to parsing the warnings text for runs annotated before that flag was
        turned on, being retried upload-only, or whose file is missing. """
    filename = annotation_run.vep_skipped_variants_filename
    if filename and os.path.exists(filename):
        skipped_count = len(parse_skipped_variants_file(filename))
        warnings_count = parse_vep_warnings(annotation_run.vep_warnings).skipped_count
        if warnings_count != skipped_count:
            logging.warning("AnnotationRun %s: --skipped_variants_file has %d rows but vep_warnings parsed "
                            "to %d skips - using the file", annotation_run.pk, skipped_count, warnings_count)
        return skipped_count, filename
    return parse_vep_warnings(annotation_run.vep_warnings).skipped_count, "vep_warnings"


def check_all_dumped_records_accounted_for(annotation_run: AnnotationRun, records_read: int):
    """ #1701: every record we handed VEP must come back either annotated or reported as skipped.

        VEP can exit 0 having written a short output VCF (Ensembl/ensembl-vep#2033 - the external
        compressor dies and Runner::run discards close()'s result), and nothing downstream notices: the
        missing records are written as vep_skipped_reason=UNKNOWN rows, which make the variants look
        annotated so they are never offered for annotation again. """
    load_vep_warnings(annotation_run)
    if annotation_run.dump_count is None:
        return  # Legacy run predating dump_count, being retried upload-only - nothing to check against

    skipped_count, source = count_vep_skipped(annotation_run)
    expected = records_read + skipped_count
    logging.info("AnnotationRun %s: %d records read + %d skipped (from %s) vs dump_count %d",
                 annotation_run.pk, records_read, skipped_count, source, annotation_run.dump_count)
    if expected != annotation_run.dump_count:
        shortfall = annotation_run.dump_count - expected
        msg = (f"AnnotationRun {annotation_run.pk}: VEP output does not account for every dumped record - "
               f"dump_count={annotation_run.dump_count}, records read={records_read}, "
               f"VEP reported {skipped_count} skipped (from {source}), shortfall {shortfall}. "
               f"The annotated VCF is truncated - retry the full run (not upload only) so it is "
               f"re-annotated from scratch.")
        raise TruncatedVEPOutputError(msg)


def handle_vep_skipped(annotation_run: AnnotationRun, bulk_inserter):
    """ set annotation_run fields 'vep_warnings' and 'vep_skipped_count'
        Create VariantAnnotation entries for any skipped variants, filling in vep_skipped_reason

        Uses bulk_inserter as that handles DB partitions """

    annotation_run.vep_skipped_count = _vep_skipped_count(annotation_run)
    if annotation_run.vep_skipped_count:
        # These will be skipped by VEP (but also we don't write out ones we know will be skipped (eg too long)
        vep_warnings = parse_vep_warnings(load_vep_warnings(annotation_run))
        incomplete_variant_ids = vep_warnings.incomplete_variant_ids
        skipped_contigs = vep_warnings.skipped_contigs

        version = annotation_run.annotation_range_lock.version
        annotation_version = version.get_any_annotation_version()
        # This pulls down any un-annotated variants (which after running VEP + inserting means were skipped)
        for v in get_variants_qs_for_annotation(annotation_version,
                                                pipeline_type=annotation_run.pipeline_type,
                                                min_variant_id=annotation_run.annotation_range_lock.min_variant_id,
                                                max_variant_id=annotation_run.annotation_range_lock.max_variant_id):
            if v.locus.contig.name in skipped_contigs:
                reason = VEPSkippedReason.UNKNOWN_CONTIG
            elif v.pk in incomplete_variant_ids:
                reason = VEPSkippedReason.INCOMPLETE
            elif v.length > settings.ANNOTATION_VEP_SV_MAX_SIZE:
                reason = VEPSkippedReason.TOO_LONG
            else:
                reason = VEPSkippedReason.UNKNOWN

            va = {"version_id": version.pk, "variant_id": v.pk,
                  "annotation_run_id": annotation_run.pk, "vep_skipped_reason": reason,
                  "predictions_num_pathogenic": 0, "predictions_num_benign": 0}
            if (annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT
                    and reason == VEPSkippedReason.TOO_LONG):
                variant_coordinate = VariantCoordinate(chrom=v.locus.contig.name,
                                                       position=v.locus.position,
                                                       ref=v.locus.ref.seq,
                                                       alt=v.alt.seq,
                                                       svlen=v.svlen)
                bulk_inserter.add_sv_gene_overlaps(v.pk, variant_coordinate, va)
            bulk_inserter.variant_annotation_list.append(va)

        # Only insert the columns we care about
        bulk_inserter.all_variant_columns = ["vep_skipped_reason"]
        bulk_inserter.finish()
