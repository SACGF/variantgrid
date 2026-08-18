import logging
import os

from django.conf import settings
from django.utils import timezone

from annotation.annotation_run_files import get_annotsv_dir
from annotation.annotsv_annotation import (
    get_annotsv_command,
    get_annotsv_command_line_version,
    get_annotsv_tsv_filename,
)
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines.base import AnnotationPipelineRunner
from annotation.vcf_files.bulk_annotsv_tsv_inserter import import_annotsv_tsv
from library.utils import execute_cmd


class AnnotSVRunner(AnnotationPipelineRunner):
    """ AnnotSV over the structural variants of a range, updating the VariantAnnotation rows the SV VEP
        run wrote. Its own pipeline rather than a stage inside that run (#720), which is what lets it be
        enabled later and backfilled, run over SVs too large for VEP, roll its bundle without a VEP
        re-annotation, and fail without taking VEP's committed data with it. """

    pipeline_type = VariantAnnotationPipelineType.ANNOTSV
    selects_unannotated = False  # takes every SV in range; depends_on guarantees VEP wrote their rows
    versioned = True  # rolling the binary or bundle re-runs the genome, without touching VEP's version
    # AnnotSV's header check requires a FORMAT column, so the dump carries a dummy sample. The genotype
    # only feeds its Samples_ID reporting.
    dump_samples = ["variantgrid"]

    def get_current_tool_version(self, genome_build) -> dict:
        return {
            "code_version": get_annotsv_command_line_version(),
            "data_version": settings.ANNOTATION_ANNOTSV_BUNDLE_VERSION,
        }

    # Deliberately no ANNOTATION_VEP_SV_MAX_SIZE filter in get_variants_qs. That cap exists because VEP
    # fills the logs with 'too long to annotate' above it; AnnotSV handles large SVs, and ranking them is
    # what it is for. Being fed the VEP dump is the reason it has never seen them.

    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        return [dump_filename, get_annotsv_tsv_filename(dump_filename, get_annotsv_dir(annotation_run))]

    def tool_finished(self, annotation_run, dump_filename) -> bool:
        # AnnotSV writes its TSV at the end of the run, so its presence is proof the run completed - no
        # separate marker needed (unlike VEP, which opens its output up front - see VEPRunner).
        return os.path.exists(get_annotsv_tsv_filename(dump_filename, get_annotsv_dir(annotation_run)))

    def record_resume_state(self, annotation_run, dump_filename):
        annotation_run.vcf_annotated_filename = get_annotsv_tsv_filename(dump_filename,
                                                                         get_annotsv_dir(annotation_run))

    def annotate(self, annotation_run, lease_heartbeat=None):
        self.check_tool_version(annotation_run)

        # #1658: per-task dump path so a reclaimed run's fresh attempt never shares one with a stalled
        # zombie. The TSV is named off the dump stem, so it is per-attempt too.
        dump_count = self.dump(annotation_run, task_token=annotation_run.task_id)
        if not dump_count:
            annotation_run.annotated_count = 0
            annotation_run.annotation_end = timezone.now()
            annotation_run.save()
            return

        annotsv_dir = get_annotsv_dir(annotation_run)
        os.makedirs(annotsv_dir, exist_ok=True)
        vcf_dump_filename = annotation_run.vcf_dump_filename
        cmd = get_annotsv_command(vcf_dump_filename, annotsv_dir,
                                  annotation_run.genome_build, annotation_run.annotation_consortium)
        annotation_run.annotation_start = timezone.now()
        annotation_run.pipeline_command = " ".join(cmd)
        annotation_run.save()

        # #1658: register AnnotSV with the lease heartbeat so a reclaim kills it rather than letting the
        # losing worker run to completion against a range the new attempt now owns.
        process_callback = lease_heartbeat.set_process if lease_heartbeat else None
        try:
            return_code, std_out, std_err = execute_cmd(
                cmd, process_callback=process_callback,
                timeout=settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS)
        finally:
            if lease_heartbeat:
                lease_heartbeat.set_process(None)  # drop the finished handle so a later tick can't kill anything

        max_output = 1_000_000
        annotation_run.pipeline_stdout = std_out[:max_output] if std_out else std_out
        annotation_run.pipeline_stderr = std_err[:max_output] if std_err else std_err
        logging.info("AnnotSV returned code: %s", return_code)

        tsv_filename = get_annotsv_tsv_filename(vcf_dump_filename, annotsv_dir)
        if return_code != 0 or not os.path.exists(tsv_filename):
            # As VEP: this is where the heartbeat's abort lands, so persist the output only while we still
            # own the run, or it goes onto the row the new attempt now holds.
            if not annotation_run.save_if_owner(annotation_run.task_id):
                logging.warning("AnnotationRun %s was reclaimed - discarding AnnotSV output from task %s",
                                annotation_run.pk, annotation_run.task_id)
            tsv_missing = "" if os.path.exists(tsv_filename) else f" (expected TSV not found: {tsv_filename})"
            raise RuntimeError(f"AnnotSV returned {return_code}{tsv_missing}")

        annotation_run.vcf_annotated_filename = tsv_filename
        annotation_run.annotation_end = timezone.now()
        annotation_run.save()

    def import_results(self, annotation_run):
        annotation_run.upload_start = timezone.now()
        annotation_run.upload_attempts += 1
        annotation_run.save()
        annotation_run.annotated_count = import_annotsv_tsv(annotation_run)
        annotation_run.upload_end = timezone.now()
        annotation_run.save()
