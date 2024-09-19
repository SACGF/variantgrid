import logging
import os
import uuid
from typing import Optional

import celery

from django.conf import settings
from django.utils import timezone

from analysis.analysis_templates import get_cohort_analysis, get_sample_analysis
from analysis.grid_export import node_grid_get_export_iterator
from analysis.models import AnalysisTemplate, SampleNode
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


def _write_cached_generated_file(cgf: CachedGeneratedFile, total_records, filename, file_iterator):
    update_size = max(1000, total_records / 100)  # 1% or every 1k records
    update_progress_iterator = update_cgf_progress_iterator(file_iterator(), cgf.pk, total_records, update_size)

    logging.info("Starting to write %s", filename)
    media_root_filename = os.path.join(settings.MEDIA_ROOT, str(uuid.uuid4()), filename)
    try:
        mk_path_for_file(media_root_filename)
        with open(media_root_filename, "w") as f:
            for line in update_progress_iterator:
                f.write(line)  # Already has newline
        cgf.filename = media_root_filename
        cgf.task_status = "SUCCESS"
        cgf.generate_end = timezone.now()
        logging.info("Wrote %s", media_root_filename)
    except Exception as e:
        logging.error("Failed to write %s: %s", media_root_filename, e)
        cgf.exception = str(e)
        cgf.task_status = "FAILURE"
    cgf.save()


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
    basename = "_".join([name_from_filename(cohort.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(cohort.genome_build)])

    request = FakeRequest(user=admin_bot())
    filename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename)
    _write_cached_generated_file(cgf, node.count, filename, file_iterator)


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
    basename = "_".join([name_from_filename(sample.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(sample.genome_build)])
    request = FakeRequest(user=admin_bot())
    filename, file_iterator = node_grid_get_export_iterator(request, node, export_type, basename=basename)
    _write_cached_generated_file(cgf, node.count, filename, file_iterator)
