import gzip
import logging
import os
import zipfile
from typing import Optional

import celery

from django.conf import settings
from django.utils import timezone

from analysis.analysis_templates import get_cohort_analysis, get_sample_analysis
from analysis.grid_export import node_grid_get_export_iterator
from analysis.models import AnalysisTemplate, SampleNode, NodeStatus
from analysis.tasks.node_update_tasks import wait_for_node
from library.django_utils import FakeRequest
from library.guardian_utils import admin_bot
from library.utils import name_from_filename, sha256sum_str, mk_path_for_file
from snpdb.models import Cohort, Sample, CachedGeneratedFile


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
            progress = i / total_records
            cgf_qs.update(progress=progress)
        yield record
    cgf_qs.update(progress=1, task_status='SUCCESS')


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
    update_progress_iterator = update_cgf_progress_iterator(file_iterator(), cgf.pk, total_records, update_size)

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


@celery.shared_task
def export_cohort_to_downloadable_file(cohort_id, export_type):
    # This should have been created in analysis.views.views_grid.cohort_grid_export
    params_hash = get_grid_downloadable_file_params_hash(cohort_id, export_type)
    cgf = CachedGeneratedFile.objects.get(generator="export_cohort_to_downloadable_file",
                                          params_hash=params_hash)

    cohort = Cohort.objects.get(pk=cohort_id)
    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
    analysis = get_cohort_analysis(cohort, analysis_template)
    node = analysis.analysisnode_set.get_subclass(output_node=True)  # Should only be 1
    if not NodeStatus.is_ready(node.status):
        wait_for_node(node.pk)  # Needs to be ready
        node = node.get_subclass()  # Easy way to reload
        if node.count is None:
            raise ValueError(f"Node {node.pk} count is None")
    _write_node_to_cached_generated_file(cgf, analysis, node, cohort.name, export_type)


@celery.shared_task
def export_sample_to_downloadable_file(sample_id, export_type):
    # This should have been created in analysis.views.views_grid.sample_grid_export
    params_hash = get_grid_downloadable_file_params_hash(sample_id, export_type)
    cgf = CachedGeneratedFile.objects.get(generator="export_sample_to_downloadable_file",
                                          params_hash=params_hash)

    sample = Sample.objects.get(pk=sample_id)
    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
    analysis = get_sample_analysis(sample, analysis_template)
    node = SampleNode.objects.get(analysis=analysis, output_node=True)  # Should only be 1
    _write_node_to_cached_generated_file(cgf, analysis, node, sample.name, export_type)
