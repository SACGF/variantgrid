from django.db import migrations

from upload.models.models_enums import UploadedFileTypes


def _add_all_fusions_file_type(apps, _schema_editor):
    """ Existing UploadSettings were created before this file type existed, so uploads of it are
        filtered out of the upload page's file list """
    UploadSettings = apps.get_model("upload", "UploadSettings")
    UploadSettingsFileType = apps.get_model("upload", "UploadSettingsFileType")

    file_type = UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS
    existing = set(UploadSettingsFileType.objects.filter(file_type=file_type)
                   .values_list("upload_settings_id", flat=True))
    records = [UploadSettingsFileType(upload_settings_id=pk, file_type=file_type)
               for pk in UploadSettings.objects.exclude(pk__in=existing).values_list("pk", flat=True)]
    UploadSettingsFileType.objects.bulk_create(records)


def _remove_all_fusions_file_type(apps, _schema_editor):
    UploadSettingsFileType = apps.get_model("upload", "UploadSettingsFileType")
    UploadSettingsFileType.objects.filter(file_type=UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("upload", "0038_alter_uploadedvcfpipelinemaxvariant_pipeline_type"),
    ]

    operations = [
        migrations.RunPython(_add_all_fusions_file_type, _remove_all_fusions_file_type),
    ]
