from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_uploaded_vcfs_needing_backfill(apps):
    """ VCFs imported before UploadedVCFPipelineMaxVariant existed have a vcf but no per-type rows,
        so the annotation-completeness check would treat them as having nothing to wait for (#1656) """
    UploadedVCF = apps.get_model("upload", "UploadedVCF")
    return UploadedVCF.objects.filter(vcf__isnull=False, pipeline_max_variants__isnull=True).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("upload", "0031_remove_uploadedvcf_max_variant_and_more"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_backfill_uploaded_vcf_pipeline_max_variant"]),
                        note="Backfill per-pipeline-type max-variant rows for existing VCFs so annotation "
                             "completeness is evaluated per pipeline type (#1656)",
                        test=_has_uploaded_vcfs_needing_backfill),
    ]
