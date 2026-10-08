from django.db import migrations


def _backfill_uploaded_combined_variant_outputs(apps, _schema_editor):
    """ CombinedVariantOutputs imported as a VCF (before #1904) have no UploadData of their own, so the
        upload page finds no data to link to - give each file that wrote a row one """
    DragenTSO500CombinedVariantOutput = apps.get_model("seqauto", "DragenTSO500CombinedVariantOutput")
    UploadedDragenTSO500CombinedVariantOutput = apps.get_model("upload", "UploadedDragenTSO500CombinedVariantOutput")

    file_upload_ids = DragenTSO500CombinedVariantOutput.objects.filter(
        file_upload__isnull=False,
        file_upload__uploadeddragentso500combinedvariantoutput__isnull=True,
    ).values_list("file_upload_id", flat=True).distinct()
    UploadedDragenTSO500CombinedVariantOutput.objects.bulk_create(
        [UploadedDragenTSO500CombinedVariantOutput(file_upload_id=file_upload_id) for file_upload_id in file_upload_ids]
    )


class Migration(migrations.Migration):

    dependencies = [
        ("upload", "0051_remove_uploadsettings"),
        ("seqauto", "0049_dragentso500combinedvariantoutput"),
    ]

    operations = [
        migrations.RunPython(_backfill_uploaded_combined_variant_outputs, migrations.RunPython.noop),
    ]
