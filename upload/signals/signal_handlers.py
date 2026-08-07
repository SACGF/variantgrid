from upload.models import UploadedVCFPendingAnnotation


def annotation_run_complete_signal_handler(variant_annotation_version, pipeline_type=None, **kwargs):
    """ Re-start any pipelines that were waiting for annotation """

    genome_build = variant_annotation_version.genome_build
    pending_qs = UploadedVCFPendingAnnotation.objects.filter(uploaded_vcf__vcf__genome_build=genome_build,
                                                             finished__isnull=True)
    if pipeline_type is not None:
        # Only VCFs containing variants of the just-completed pipeline type could have been unblocked by it
        pending_qs = pending_qs.filter(uploaded_vcf__pipeline_max_variants__pipeline_type=pipeline_type)
    for pending_annotation in pending_qs.distinct():
        # is_fully_annotated() makes the final per-type call
        pending_annotation.attempt_schedule_annotation_stage_steps(variant_annotation_version)
