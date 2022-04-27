from annotation.annotation_versions import get_lowest_unannotated_variant_id
from upload.models import UploadedVCFPendingAnnotation, UploadSettings


def annotation_run_complete_signal_handler(variant_annotation_version, **kwargs):
    """ Re-start any pipelines that were waiting for annotation """

    lowest_unannotated_variant = get_lowest_unannotated_variant_id(variant_annotation_version)
    genome_build = variant_annotation_version.genome_build
    for pending_annotation in UploadedVCFPendingAnnotation.objects.filter(uploaded_vcf__max_variant_id__lt=lowest_unannotated_variant,
                                                                          uploaded_vcf__vcf__genome_build=genome_build,
                                                                          finished__isnull=True):
        pending_annotation.attempt_schedule_annotation_stage_steps(lowest_unannotated_variant)

    # TODO: Check for any hung pipelines where there's no other annotation runs?


def upload_settings_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        UploadSettings.create_default_visible_file_types()

