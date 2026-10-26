from eventlog.models import api_request_note
from upload.models import BackendVCF, UploadedVCFPendingAnnotation
from upload.vcf.vcf_import import link_samples_and_vcfs_to_sequencing


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


def current_sample_sheet_data_moved_handler(sample_sheet, relink_samples, **kwargs):
    """ Relink the samples of the sheet's imported joint called VCFs - all of them when relink_samples
        (members moved, or a sequencing sample is unlinked), otherwise those that need it """
    relinked_vcfs = []
    for joint_called_vcf in sample_sheet.get_joint_called_vcfs():
        try:
            backend_vcf = joint_called_vcf.backendvcf
        except BackendVCF.DoesNotExist:
            continue

        if relink_samples or joint_called_vcf.needs_to_be_linked():
            try:
                link_samples_and_vcfs_to_sequencing(backend_vcf, replace_existing=True)
                relinked_vcfs.append(str(backend_vcf.vcf))
            except ValueError as e:
                api_request_note(f"Could not relink samples of {backend_vcf.vcf}: {e}")
    if relinked_vcfs:
        api_request_note(f"Relinked VCF samples to sequencing samples: {', '.join(relinked_vcfs)}")
