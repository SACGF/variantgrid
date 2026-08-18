"""
The one place an AnnotationRun's on-disk output is removed (#1670).

Every artifact the pipeline writes - the variant dump, the VEP annotated VCF, the conservation sidecar
and the AnnotSV directory - is reclaimed here and nowhere else. Callers do not remove files themselves;
they announce what happened and this module decides, so the "which files does a run own" logic has a
single owner and the pipeline modules stay out of the filesystem.

Three ways in, one exit:

  annotation_run_complete_signal   the run's rows are committed, so its output has done its job.
                                   Gated on ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS.
  annotation_run_discarded_signal  an attempt's output is being thrown away (a losing attempt
                                   unwinding, or a stalled run scrubbed back to CREATED). Always
                                   removed - it is unreachable, and a re-dump would collide with it.
  post_delete                      the row is going away, and it was the only thing that could name
                                   these files. Using the signal rather than overriding
                                   AnnotationRun.delete() also covers queryset deletes (registering a
                                   receiver disables Django's fast-delete path), which an override
                                   would silently miss.

Deliberately absent: a failure path. A run that errored keeps everything on disk for investigation,
which is the behaviour #1670 asked for - so it is achieved by there being no signal to send.
"""
import logging
import os
import shutil

from django.conf import settings

from annotation.annotation_run_files import get_annotsv_dir
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines import get_runner
from annotation.sv_conservation import conservation_sidecar_filename
from annotation.vep_annotation import get_vep_skipped_variants_filename


def get_all_run_output_paths(annotation_run) -> list[str]:
    """ Every path this run may own, from the derived dump stem *and* the persisted filename fields.

        The two sources cover each other. Derivation (#1660) names files a run was interrupted before it
        could record - vcf_annotated_filename only reaches the DB at the final save in the runner's
        annotate(), while the dump stem is persisted before the tool starts. The persisted fields in turn
        name files whose path diverges from the derivation: an external run (#1568) imported without a
        local dump has no dump stem at all, so only the field names its annotated VCF. """
    paths = []
    if dump_filename := annotation_run.vcf_dump_filename:
        paths += get_runner(annotation_run.pipeline_type).get_output_paths(annotation_run, dump_filename)
    if annotated_filename := annotation_run.vcf_annotated_filename:
        paths += [annotated_filename, conservation_sidecar_filename(annotated_filename),
                  get_vep_skipped_variants_filename(annotated_filename)]
    if skipped_variants_filename := annotation_run.vep_skipped_variants_filename:
        paths.append(skipped_variants_filename)
    return paths


def remove_annotation_run_output(annotation_run, remove_annotsv_dir: bool = True):
    """ Best-effort removal of everything this run owns on disk.

        remove_annotsv_dir=False is for a caller holding only its *own* attempt's files: the AnnotSV
        directory is keyed on the run and shared between attempts (get_annotsv_dir), so a discarded
        attempt deleting it would take the winning attempt's output with it. A run that is finished or
        being deleted has no such rival, and clearing the whole directory is the only way the non-TSV
        files AnnotSV writes alongside its output are ever reclaimed. """
    for path in dict.fromkeys(get_all_run_output_paths(annotation_run)):  # de-dupe, keeping order
        if path and os.path.exists(path):
            try:
                os.remove(path)
                logging.info("Removed AnnotationRun %s output file: %s", annotation_run.pk, path)
            except OSError:
                logging.exception("Failed removing AnnotationRun %s output file: %s",
                                  annotation_run.pk, path)
    if remove_annotsv_dir and annotation_run.pipeline_type == VariantAnnotationPipelineType.ANNOTSV:
        annotsv_dir = get_annotsv_dir(annotation_run)
        if os.path.isdir(annotsv_dir):
            logging.info("Removing AnnotationRun %s AnnotSV dir: %s", annotation_run.pk, annotsv_dir)
            shutil.rmtree(annotsv_dir, ignore_errors=True)


def annotation_run_complete_cleanup_handler(sender, annotation_run=None, **kwargs):  # pylint: disable=unused-argument
    """ The run imported successfully - its output has served its purpose, so reclaim the disk.

        annotation_run is optional because this signal predates the cleanup and other senders may not
        carry it. The filename fields are deliberately left pointing at the removed paths as a record of
        what the run produced: NULLing vcf_annotated_filename would move an external run's get_status()
        off FINISHED (see AnnotationRun.get_status). """
    if annotation_run is None:
        return
    if not settings.ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS:
        return
    remove_annotation_run_output(annotation_run)


def annotation_run_discarded_cleanup_handler(sender, annotation_run, **kwargs):  # pylint: disable=unused-argument
    """ This attempt's output is being thrown away - remove it, but leave the shared AnnotSV directory. """
    remove_annotation_run_output(annotation_run, remove_annotsv_dir=False)


def annotation_run_post_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    """ The row is gone and nothing else can name these files - there is no sweep of
        ANNOTATION_VCF_DUMP_DIR to catch them later. """
    remove_annotation_run_output(instance)
