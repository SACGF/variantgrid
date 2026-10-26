"""
Refreshes ResolvedVariantInfo.p_hgvs once VEP has annotated its variant (#2421) - classifications are matched straight
after variant insert, before VEP has run, so they start on the computed / ClinGen / evidence fallbacks.

Entry point: annotation_run_complete_p_hgvs_handler (annotation_run_complete_signal), which queues
recalc_p_hgvs_for_annotated_variants_task so the annotation lane never waits on it.
"""
from typing import Optional

from django.db import transaction
from django.dispatch import receiver

from annotation.models.models_annotation_run import AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.pipelines import vep_pipeline_types
from annotation.signals.manual_signals import annotation_run_complete_signal
from classification.tasks.resolved_variant_info_p_hgvs_task import (
    recalc_p_hgvs_for_annotated_variants_task,
)


@receiver(annotation_run_complete_signal)
def annotation_run_complete_p_hgvs_handler(sender, variant_annotation_version: VariantAnnotationVersion,
                                           pipeline_type=None, annotation_run: Optional[AnnotationRun] = None,
                                           **kwargs):  # pylint: disable=unused-argument
    if annotation_run is None or pipeline_type not in vep_pipeline_types():
        return
    genome_build = variant_annotation_version.genome_build
    # p.HGVS reads hgvs_p from the latest version only - a new version being annotated is picked up by
    # the recalc_resolved_variant_info_p_hgvs command once it becomes the latest
    if variant_annotation_version != VariantAnnotationVersion.latest(genome_build):
        return
    range_lock = annotation_run.annotation_range_lock
    task = recalc_p_hgvs_for_annotated_variants_task.si(genome_build.name, range_lock.min_variant_id,
                                                        range_lock.max_variant_id)
    transaction.on_commit(task.apply_async)
