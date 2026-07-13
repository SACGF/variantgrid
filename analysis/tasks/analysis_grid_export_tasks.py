import gzip
import logging
import os
import zipfile
from typing import Optional

import celery
from celery.exceptions import Retry
from django.conf import settings
from django.utils import timezone

from analysis.analysis_templates import get_cohort_analysis, get_sample_analysis
from analysis.grid_export import node_grid_get_export_iterator
from analysis.models import AnalysisTemplate, NodeStatus, CohortNode, SampleNode
from library.constants import MINUTE_SECS
from library.django_utils import FakeRequest
from library.guardian_utils import admin_bot
from library.log_utils import log_traceback
from library.utils import name_from_filename, sha256sum_str, mk_path_for_file
from snpdb.models import Cohort, Sample, CachedGeneratedFile

# How long (secs) to wait between checks that an export's output node has finished loading.
# We re-queue the export task (Celery retry) between checks rather than sleeping, so the worker
# is freed and we don't deadlock when workers are busy loading the parent nodes we're waiting on.
NODE_WAIT_TIME_BETWEEN_CHECKS = [5, 5, 10, 10, 30, 30, 60, MINUTE_SECS * 2]


def get_annotated_download_files_cgf(generator, pk) -> dict[str, Optional[CachedGeneratedFile]]:
    annotated_download_files = {}
    try:
        AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        params_hash_vcf = get_grid_downloadable_file_params_hash(pk, "vcf")
        cgf_vcf = CachedGeneratedFile.objects.filter(generator=generator,
                                                     params_hash=params_hash_vcf).first()
        params_hash_csv = get_grid_downloadable_file_params_hash(pk, "csv")
        cgf_csv = CachedGeneratedFile.objects.filter(generator=generator,
                                                     params_hash=params_hash_csv).first()

        annotated_download_files = {"vcf": cgf_vcf, "csv": cgf_csv}
    except ValueError:
        pass

    return annotated_download_files


def get_grid_downloadable_file_params_hash(pk, export_type):
    return sha256sum_str(f"{pk}-{export_type}")


def update_cgf_progress_iterator(iterator, cgf_id, total_records, update_size):
    update_size = int(update_size)  # make sure int so modulus below will hit
    cgf_qs = CachedGeneratedFile.objects.filter(id=cgf_id)
    cgf_qs.update(progress=0)

    for i, record in enumerate(iterator):
        if i % update_size == 0:
            progress = i / total_records if total_records else 0
            cgf_qs.update(progress=progress)
        yield record
    cgf_qs.update(progress=1)


def _write_node_to_cached_generated_file(cgf, analysis, node, name, export_type):
    basename = "_".join([name_from_filename(name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(analysis.genome_build)])
    request = FakeRequest(user=admin_bot())
    basename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename)
    open_func = open
    if export_type == 'vcf':
        open_func = gzip.open
        basename += ".gz"

    total_records = node.count
    update_size = max(1000, total_records / 100)  # 1% or every 1k records
    update_progress_iterator = update_cgf_progress_iterator(file_iterator, cgf.pk, total_records, update_size)

    media_root_filename = os.path.join(settings.GENERATED_DIR, cgf.generator, str(cgf.pk), basename)
    logging.info("Starting to write %s", media_root_filename)
    try:
        mk_path_for_file(media_root_filename)
        with open_func(media_root_filename, "wt") as f:
            for line in update_progress_iterator:
                f.write(line)  # Already has newline

        if export_type == 'csv':
            original_filename = media_root_filename
            zip_file_path = media_root_filename + ".zip"
            with zipfile.ZipFile(zip_file_path, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(original_filename, arcname=os.path.basename(original_filename))
            os.unlink(original_filename)
            media_root_filename = zip_file_path
        cgf.filename = media_root_filename
        cgf.task_status = "SUCCESS"
        cgf.generate_end = timezone.now()
        logging.info("Wrote %s", media_root_filename)
        # Write CSVs to Zip (requires the file to be there already)
    except Exception as e:
        logging.error("Failed to write %s: %s", media_root_filename, e)
        cgf.exception = str(e)
        cgf.task_status = "FAILURE"
    cgf.save()

def _wait_for_output_node(self, node):
    """ Ensure the export's output node is ready.

        `self` is the bound export task. If the node is still loading we re-queue the export task
        (self.retry) rather than blocking the worker - each retry re-runs the export from the top,
        re-fetching the node. Raises once the node is ready but has no count, or once we've waited
        the full NODE_WAIT_TIME_BETWEEN_CHECKS schedule without it becoming ready. """
    if not NodeStatus.is_ready(node.status):
        retry_index = self.request.retries
        if retry_index < len(NODE_WAIT_TIME_BETWEEN_CHECKS):
            countdown = NODE_WAIT_TIME_BETWEEN_CHECKS[retry_index]
            logging.info("Output node %s not ready (status=%s) - retrying export in %d secs (attempt %d)",
                         node.pk, node.get_status_display(), countdown, retry_index + 1)
            raise self.retry(countdown=countdown, max_retries=len(NODE_WAIT_TIME_BETWEEN_CHECKS))
        raise ValueError(f"Node {node}/{node.version} status={node.get_status_display()} "
                         "still not ready after waiting - aborting export")
    node = node.get_subclass()  # Easy way to reload
    if node.count is None:
        raise ValueError(f"Node {node}/{node.version} - {node.status} count is None")
    return node

@celery.shared_task(bind=True)
def export_cohort_to_downloadable_file(self, cohort_id, export_type):
    try:
        # This should have been created in analysis.views.views_grid.cohort_grid_export
        params_hash = get_grid_downloadable_file_params_hash(cohort_id, export_type)
        cgf = CachedGeneratedFile.objects.get(generator="export_cohort_to_downloadable_file",
                                              params_hash=params_hash)

        cohort = Cohort.objects.get(pk=cohort_id)
        analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        analysis = get_cohort_analysis(cohort, analysis_template)
        node = CohortNode.objects.get_subclass(analysis=analysis, output_node=True)  # Should only be 1
        node = _wait_for_output_node(self, node)
        _write_node_to_cached_generated_file(cgf, analysis, node, cohort.name, export_type)
    except Retry:
        raise  # Output node not ready yet - export task re-queued, not an error
    except Exception:
        log_traceback()
        raise

@celery.shared_task(bind=True)
def export_sample_to_downloadable_file(self, sample_id, export_type):
    try:
        # This should have been created in analysis.views.views_grid.sample_grid_export
        params_hash = get_grid_downloadable_file_params_hash(sample_id, export_type)
        cgf = CachedGeneratedFile.objects.get(generator="export_sample_to_downloadable_file",
                                              params_hash=params_hash)

        sample = Sample.objects.get(pk=sample_id)
        analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
        analysis = get_sample_analysis(sample, analysis_template)
        node = SampleNode.objects.get_subclass(analysis=analysis, output_node=True)  # Should only be 1
        node = _wait_for_output_node(self, node)
        _write_node_to_cached_generated_file(cgf, analysis, node, sample.name, export_type)
    except Retry:
        raise  # Output node not ready yet - export task re-queued, not an error
    except Exception:
        log_traceback()
        raise
