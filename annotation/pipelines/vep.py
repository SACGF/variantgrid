import logging
import os

from django.conf import settings
from django.db.models.functions.math import Abs
from django.db.models.query_utils import Q
from django.utils import timezone

from annotation.annotation_run_files import get_annotated_filename
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines.base import AnnotationPipelineRunner
from annotation.sv_conservation import (
    conservation_sidecar_filename,
    get_sv_conservation_tracks,
    score_sv_vcf,
    write_conservation_sidecar,
)
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import (
    get_vep_command,
    get_vep_skipped_variants_filename,
    vep_check_command_line_version_match,
)
from library.log_utils import log_traceback
from library.utils import execute_cmd


class VEPRunner(AnnotationPipelineRunner):
    """ Runs VEP over a class of variant. One instance per pipeline type - the class of variant is the
        only thing that differs between them, and get_vep_command already takes it. """

    def __init__(self, pipeline_type: VariantAnnotationPipelineType):
        self.pipeline_type = pipeline_type

    @property
    def is_structural_variant(self) -> bool:
        return self.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT

    def get_variants_qs(self, annotation_run):
        qs = super().get_variants_qs(annotation_run)
        if self.is_structural_variant and settings.ANNOTATION_VEP_SV_MAX_SIZE:
            # VEP will skip variants above a certain size and fill up the logs with 'too long to annotate'
            # So just skip these. I don't think it makes much difference in memory usage
            q_not_too_long = Q(svlen__isnull=True) | Q(abs_svlen__lte=settings.ANNOTATION_VEP_SV_MAX_SIZE)
            qs = qs.annotate(abs_svlen=Abs("svlen")).filter(q_not_too_long)
        return qs

    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        annotated_filename = get_annotated_filename(annotation_run, dump_filename)
        return [
            dump_filename,
            annotated_filename,
            conservation_sidecar_filename(annotated_filename),
            get_vep_skipped_variants_filename(annotated_filename),
        ]

    def tool_finished(self, annotation_run, dump_filename) -> bool:
        # #1710: VEP creates the annotated VCF at startup and writes to it progressively, so its presence
        # says nothing about how far the run got - a run reclaimed mid-VEP used to resume upload-only off
        # a part-written file, and the #1701 checks in the import lane then failed it. Runner::run closes
        # the output handle and only then calls finish(), which writes the skipped-variants list (even
        # when nothing was skipped), so that file appearing is the proof VEP reached the end. It is named
        # off the per-attempt dump stem, so it can only have come from this attempt.
        annotated_filename = get_annotated_filename(annotation_run, dump_filename)
        if not os.path.exists(get_vep_skipped_variants_filename(annotated_filename)):
            return False
        return os.path.exists(annotated_filename)

    def record_resume_state(self, annotation_run, dump_filename):
        annotated_filename = get_annotated_filename(annotation_run, dump_filename)
        annotation_run.vcf_annotated_filename = annotated_filename
        # #1701: record VEP's skipped-variants list so the upload-only relaunch checks records-in vs
        # records-out against VEP's own list rather than the warnings-text fallback.
        annotation_run.vep_skipped_variants_filename = get_vep_skipped_variants_filename(annotated_filename)

    def annotate(self, annotation_run, lease_heartbeat=None, vep_version_check=True):
        if vep_version_check:
            # Do a check before we annotate
            vep_check_command_line_version_match(annotation_run.variant_annotation_version)

        # #1658: tag the local dump (and, via the derived annotated name, the annotated VCF) with the
        # executing task_id so a reclaimed run's fresh attempt never shares a path with a stalled zombie.
        vcf_dump_count = self.dump(annotation_run, task_token=annotation_run.task_id)
        vcf_dump_filename = annotation_run.vcf_dump_filename

        genome_build = annotation_run.genome_build
        annotation_consortium = annotation_run.annotation_consortium

        logging.info("Annotating %d variants", vcf_dump_count)
        if vcf_dump_count:
            vcf_annotated_filename = get_annotated_filename(annotation_run, vcf_dump_filename)

            cmd = get_vep_command(vcf_dump_filename, vcf_annotated_filename, genome_build, annotation_consortium,
                                  annotation_run.pipeline_type,
                                  variant_annotation_version=annotation_run.variant_annotation_version)
            annotation_run.annotation_start = timezone.now()
            annotation_run.pipeline_command = " ".join(cmd)
            annotation_run.save()
            # #1658: register VEP with the lease heartbeat so a reclaim (lost lease) kills it mid-run instead
            # of letting the losing worker annotate to completion and double-write the reassigned run's state.
            process_callback = lease_heartbeat.set_process if lease_heartbeat else None
            try:
                return_code, std_out, std_err = execute_cmd(cmd, process_callback=process_callback)
            finally:
                if lease_heartbeat:
                    lease_heartbeat.set_process(None)  # drop the finished handle so a later tick can't kill anything
            # VEP can produce enormous output (>1GB) for large batches - PostgreSQL has a 1GB field limit
            max_output = 1_000_000
            annotation_run.pipeline_stdout = std_out[:max_output] if std_out else std_out
            annotation_run.pipeline_stderr = std_err[:max_output] if std_err else std_err
            logging.info(f"VEP returned code: {return_code}")

            if return_code != 0:
                # #1658: this is the path the heartbeat's abort itself lands on, so it is exactly where we may
                # no longer own the run - persist the VEP output only while we do, or it goes onto the row the
                # new attempt now holds. Guarded on our own task_id (None matches an unleased run, which is how
                # the annotation_external command calls this).
                if not annotation_run.save_if_owner(annotation_run.task_id):  # save stdout/stderr
                    logging.warning("AnnotationRun %s was reclaimed - discarding VEP output from task %s",
                                    annotation_run.pk, annotation_run.task_id)
                raise RuntimeError(f"VEP returned {return_code}")

            annotation_run.vcf_annotated_filename = vcf_annotated_filename
            # #1701: VEP writes this in Runner::finish(), after the output handle is closed - so it is present
            # and complete even for the truncated-output failure this check exists to catch. Recorded only when
            # written, so the import lane falls back to the warnings text for anything that didn't produce one.
            skipped_variants_filename = get_vep_skipped_variants_filename(vcf_annotated_filename)
            if os.path.exists(skipped_variants_filename):
                annotation_run.vep_skipped_variants_filename = skipped_variants_filename
            annotation_run.annotation_end = timezone.now()

            if settings.ANNOTATION_VEP_SV_CONSERVATION_PYBIGWIG_ENABLED and self.is_structural_variant:
                self._write_conservation_sidecar(annotation_run, vcf_dump_filename, vcf_annotated_filename)
        else:
            # Now we have standard/CNV type pipelines, it's possible some can be empty
            annotation_run.annotated_count = 0
            annotation_run.annotation_end = timezone.now()

        annotation_run.save()

    def _write_conservation_sidecar(self, annotation_run, vcf_dump_filename, vcf_annotated_filename):
        """ Conservation (phastCons/phyloP) _max columns for SVs - computed with pyBigWig instead of the
            4 conservation VEP --custom bigWig overlaps, whose O(SV-span) cost makes large SVs never finish
            (#1657). The values are written to a sidecar TSV next to the annotated VCF; the import lane
            (BulkVEPVCFAnnotationInserter) merges them into the same _max columns. Best-effort: a scoring
            failure logs and leaves the columns null rather than failing the run. """
        genome_build = annotation_run.genome_build
        try:
            tracks = get_sv_conservation_tracks(genome_build)
            results = score_sv_vcf(vcf_dump_filename, genome_build)
            sidecar = conservation_sidecar_filename(vcf_annotated_filename)
            write_conservation_sidecar(sidecar, results, tracks)
            logging.info("Wrote SV conservation for %d variants to %s", len(results), sidecar)
        except Exception:
            log_traceback()
            logging.warning("SV conservation (pyBigWig) stage failed for AnnotationRun %s",
                            annotation_run.pk)

    def import_results(self, annotation_run):
        import_vcf_annotations(annotation_run)
