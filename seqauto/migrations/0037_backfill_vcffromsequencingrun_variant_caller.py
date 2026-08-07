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
        # upload.0030 renames BackendVCF.vcf_file -> single_sample_vcf. This migration only
        # depends on upload.0029 (pre-rename), so depending on the topological order the field
        # is reached under either name. Both are null=True OneToOnes, so getattr returns None
        # when the relation (or the field itself) is absent - handle whichever name applies.
        backing = (getattr(backend_vcf, "single_sample_vcf", None)
                   or getattr(backend_vcf, "vcf_file", None)
                   or backend_vcf.joint_called_vcf)
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
