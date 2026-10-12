import gzip
import json
import logging
import os
import zipfile
from collections import defaultdict
from datetime import timedelta
from typing import Optional

import celery
from celery.exceptions import Retry
from django.conf import settings
from django.contrib.auth.models import User
from django.utils import timezone

from analysis.analysis_templates import get_cohort_analysis, get_sample_analysis
from analysis.grid_export import get_node_export_basename, node_grid_get_export_iterator
from analysis.models import AnalysisTemplate, CohortNode, NodeStatus, SampleNode, VariantTag
from analysis.models.nodes.analysis_node import node_query_planner_settings
from analysis.views.analysis_permissions import get_node_subclass_or_non_fatal_exception
from genes.models import CanonicalTranscriptCollection
from library.constants import MINUTE_SECS
from library.django_utils import FakeRequest
from library.guardian_utils import admin_bot
from library.log_utils import log_traceback, report_exc_info
from library.utils import mk_path_for_file, name_from_filename, sha256sum_str
from snpdb.models import CachedGeneratedFile, Cohort, Sample
from snpdb.utils import get_tag_sort_order_by_tag

# How long (secs) to wait between checks that an export's output node has finished loading.
# We re-queue the export task (Celery retry) between checks rather than sleeping, so the worker
# is freed and we don't deadlock when workers are busy loading the parent nodes we're waiting on.
NODE_WAIT_TIME_BETWEEN_CHECKS = [5, 5, 10, 10, 30, 30, 60, MINUTE_SECS * 2]

NODE_EXPORT_GENERATOR = "export_node_to_downloadable_file"


def _get_usable_cgf(generator, pk, export_type) -> Optional[CachedGeneratedFile]:
    """ A cached file that needs regenerating is no better than never having generated it - return None so
        the page offers the generate link, and get_or_create_and_launch drops the row when they click.
        A failed row is returned so the page shows the failure rather than a link that would fail again """
    params_hash = get_grid_downloadable_file_params_hash(pk, export_type)
    cgf = CachedGeneratedFile.objects.filter(generator=generator, params_hash=params_hash).first()
    if cgf and cgf.needs_regenerating:
        cgf = None
    return cgf


def get_annotated_download_files_cgf(generator, pk) -> dict[str, Optional[CachedGeneratedFile]]:
    annotated_download_files = {}
    try:
        AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        annotated_download_files = {
            "vcf": _get_usable_cgf(generator, pk, "vcf"),
            "csv": _get_usable_cgf(generator, pk, "csv"),
        }
    except ValueError:
        pass

    return annotated_download_files


def get_grid_downloadable_file_params_hash(pk, export_type):
    return sha256sum_str(f"{pk}-{export_type}")


def get_node_grid_downloadable_file_params_hash(node_id, node_version, user_id, export_type,
                                                canonical_transcript_collection_id, grid_params: dict) -> str:
    """ Node grid content is user-dependent (tags from other analyses) and filter-dependent, and the node
        version invalidates the cache when the analysis is edited - so all of it goes into the hash """
    params = json.dumps([node_id, node_version, user_id, export_type, canonical_transcript_collection_id,
                         grid_params], sort_keys=True)
    return sha256sum_str(params)


def update_cgf_progress_iterator(iterator, cgf_id, total_records, update_size):
    """ Wraps the export's rows iterator - the file iterator yields a chunk of rows at a time, so
        counting its yields would under-report progress by the chunk size. Each update bumps modified,
        the heartbeat that keeps a long export from looking stalled (@see CachedGeneratedFile.generation_stalled) """
    update_size = int(update_size)  # make sure int so modulus below will hit
    cgf_qs = CachedGeneratedFile.objects.filter(id=cgf_id)
    cgf_qs.update(progress=0, modified=timezone.now())

    for i, record in enumerate(iterator):
        if i % update_size == 0:
            # Capped as the cached node count can be behind the rows actually exported
            progress = min(i / total_records, 1) if total_records else 0
            cgf_qs.update(progress=progress, modified=timezone.now())
        yield record
    cgf_qs.update(progress=1, modified=timezone.now())


def _get_cgf_for_task(task, generator, params_hash) -> Optional[CachedGeneratedFile]:
    """ The row this run writes to, with its heartbeat bumped - or None if the row has been relaunched under
        another task since this one was queued (@see CachedGeneratedFile.generation_stalled), so two runs
        don't write the same file """
    # This should have been created by CachedGeneratedFile.get_or_create_and_launch
    cgf = CachedGeneratedFile.objects.get(generator=generator, params_hash=params_hash)
    if cgf.task_id and cgf.task_id != task.request.id:
        logging.info("CachedGeneratedFile %s now belongs to task %s - skipping superseded task %s",
                     cgf.pk, cgf.task_id, task.request.id)
        return None
    CachedGeneratedFile.objects.filter(pk=cgf.pk).update(modified=timezone.now())
    return cgf


def _get_annotated_basename(analysis, name: str) -> str:
    return "_".join([name_from_filename(name), "annotated", f"v{analysis.annotation_version.pk}",
                     str(analysis.genome_build)])


def _write_node_to_cached_generated_file(cgf, request, node, basename, export_type, zip_csv=True, **export_kwargs):
    """ zip_csv: whole cohort/sample exports run to millions of rows so are zipped; a node export is what the
        user filtered the grid down to, and they open it straight away, so it stays a plain .csv """
    total_records = node.count
    update_size = max(1000, total_records / 100)  # 1% or every 1k records

    def row_wrapper(rows):
        return update_cgf_progress_iterator(rows, cgf.pk, total_records, update_size)

    basename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename,
                                                            row_wrapper=row_wrapper, **export_kwargs)
    open_func = open
    if export_type == 'vcf':
        open_func = gzip.open
        basename += ".gz"

    media_root_filename = os.path.join(settings.GENERATED_DIR, cgf.generator, str(cgf.pk), basename)
    logging.info("Starting to write %s", media_root_filename)
    try:
        mk_path_for_file(media_root_filename)
        with open_func(media_root_filename, "wt") as f:
            for chunk in file_iterator:
                f.write(chunk)  # Already has newline

        if export_type == 'csv' and zip_csv:
            original_filename = media_root_filename
            zip_file_path = media_root_filename + ".zip"
            with zipfile.ZipFile(zip_file_path, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(original_filename, arcname=os.path.basename(original_filename))
            os.unlink(original_filename)
            media_root_filename = zip_file_path
        cgf.filename = media_root_filename
        cgf.task_status = "SUCCESS"
        cgf.progress = 1  # row_wrapper updated the DB directly, so our copy is still at its starting value
        cgf.generate_end = timezone.now()
        logging.info("Wrote %s", media_root_filename)  # lgtm[py/clear-text-logging-sensitive-data]
    except Exception as e:
        # The row stays failed until the cause is fixed (@see CachedGeneratedFile.needs_regenerating), so
        # admins hear about it now rather than from the user
        logging.error("Failed to write %s: %s", media_root_filename, e)  # lgtm[py/clear-text-logging-sensitive-data]
        report_exc_info()
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
        params_hash = get_grid_downloadable_file_params_hash(cohort_id, export_type)
        cgf = _get_cgf_for_task(self, "export_cohort_to_downloadable_file", params_hash)
        if cgf is None:
            return

        cohort = Cohort.objects.get(pk=cohort_id)
        analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        analysis = get_cohort_analysis(cohort, analysis_template)
        node = CohortNode.objects.get_subclass(analysis=analysis, output_node=True)  # Should only be 1
        node = _wait_for_output_node(self, node)
        request = FakeRequest(user=admin_bot())
        _write_node_to_cached_generated_file(cgf, request, node, _get_annotated_basename(analysis, cohort.name),
                                             export_type)
    except Retry:
        raise  # Output node not ready yet - export task re-queued, not an error
    except Exception:
        log_traceback()
        raise

@celery.shared_task(bind=True)
def export_sample_to_downloadable_file(self, sample_id, export_type):
    try:
        params_hash = get_grid_downloadable_file_params_hash(sample_id, export_type)
        cgf = _get_cgf_for_task(self, "export_sample_to_downloadable_file", params_hash)
        if cgf is None:
            return

        sample = Sample.objects.get(pk=sample_id)
        analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
        analysis = get_sample_analysis(sample, analysis_template)
        node = SampleNode.objects.get_subclass(analysis=analysis, output_node=True)  # Should only be 1
        node = _wait_for_output_node(self, node)
        request = FakeRequest(user=admin_bot())
        _write_node_to_cached_generated_file(cgf, request, node, _get_annotated_basename(analysis, sample.name),
                                             export_type)
    except Retry:
        raise  # Output node not ready yet - export task re-queued, not an error
    except Exception:
        log_traceback()
        raise


@celery.shared_task(bind=True)
def export_node_to_downloadable_file(self, node_id, node_version, user_id, export_type,
                                     canonical_transcript_collection_id, grid_params):
    """ Analysis node CSV/VCF export (issue #1257) - runs under Celery so a big node isn't bound by the
        gunicorn request timeout. grid_params are the whitelisted grid GET params the browser sent, so
        the export sees the same filters as the grid the user was looking at. """
    try:
        params_hash = get_node_grid_downloadable_file_params_hash(node_id, node_version, user_id, export_type,
                                                                 canonical_transcript_collection_id, grid_params)
        cgf = _get_cgf_for_task(self, NODE_EXPORT_GENERATOR, params_hash)
        if cgf is None:
            return

        user = User.objects.get(pk=user_id)
        node = get_node_subclass_or_non_fatal_exception(user, node_id, version=node_version)
        node = _wait_for_output_node(self, node)

        canonical_transcript_collection = None
        if canonical_transcript_collection_id:
            canonical_transcript_collection = CanonicalTranscriptCollection.objects.get(
                pk=canonical_transcript_collection_id)

        sort_order_by_tag = get_tag_sort_order_by_tag(user)
        # A variant carries a tag once per sample (@see VariantTag.Meta) - the export lists it once
        variant_tag_ids = defaultdict(set)
        for variant_id, tag_id in VariantTag.objects.filter(analysis=node.analysis).values_list("variant_id", "tag_id"):
            variant_tag_ids[variant_id].add(tag_id)
        variant_tags_dict = {
            variant_id: ", ".join(sorted(tag_ids, key=lambda t: (sort_order_by_tag.get(t, 0), t)))
            for variant_id, tag_ids in variant_tag_ids.items()
        }

        request = FakeRequest(user=user)
        request.GET = grid_params
        with node_query_planner_settings():
            _write_node_to_cached_generated_file(cgf, request, node, get_node_export_basename(node), export_type,
                                                 zip_csv=False,
                                                 canonical_transcript_collection=canonical_transcript_collection,
                                                 variant_tags_dict=variant_tags_dict)
    except Retry:
        raise  # Output node not ready yet - export task re-queued, not an error
    except Exception:
        log_traceback()
        raise


@celery.shared_task
def clear_old_node_export_cached_generated_files():
    """ Node exports are cached per (node, version, user, filters, export type) so build up far faster than
        the cohort/sample ones - the pre_delete signal unlinks each file as its row goes """
    cutoff = timezone.now() - timedelta(days=settings.ANALYSIS_NODE_EXPORT_CACHE_DAYS)
    CachedGeneratedFile.objects.filter(generator=NODE_EXPORT_GENERATOR, modified__lt=cutoff).delete()
