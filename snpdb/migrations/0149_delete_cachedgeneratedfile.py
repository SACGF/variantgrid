# Generated by Django 4.2.11 on 2024-09-19 04:07

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0148_one_off_move_cached_generated_files_to_uuid"),
    ]

    operations = [
        migrations.DeleteModel(
            name="CachedGeneratedFile",
        ),
    ]
