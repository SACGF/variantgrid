from django.db import migrations


def backfill_variant_caller(apps, schema_editor):
    VCFFromSequencingRun = apps.get_model("seqauto", "VCFFromSequencingRun")
    BackendVCF = apps.get_model("upload", "BackendVCF")

    rows = VCFFromSequencingRun.objects.filter(variant_caller__isnull=True)
    for row in rows.iterator():
        try:
            backend_vcf = BackendVCF.objects.get(uploaded_vcf__vcf=row.vcf)
        except BackendVCF.DoesNotExist:
            continue
        backing = backend_vcf.vcf_file or backend_vcf.joint_called_vcf
        if backing and backing.variant_caller_id is not None:
            row.variant_caller_id = backing.variant_caller_id
            row.save(update_fields=["variant_caller"])


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0036_rename_samplesheetcombinedvcffile_jointcalledvcf_and_more"),
        ("upload", "0029_rename_combo_vcf_backendvcf_joint_called_vcf"),
    ]

    operations = [
        migrations.RunPython(backfill_variant_caller, migrations.RunPython.noop),
    ]
