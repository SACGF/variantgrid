from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("upload", "0050_modifiedimportedvariant_old_value_indexes"),
    ]

    operations = [
        migrations.DeleteModel(
            name="UploadSettingsFileType",
        ),
        migrations.DeleteModel(
            name="UploadSettings",
        ),
    ]
