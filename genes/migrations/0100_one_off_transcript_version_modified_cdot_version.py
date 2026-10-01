from django.db import migrations
from django.db.models import Max

BATCH_SIZE = 100_000


def _move_data_cdot_to_modified_cdot_version(apps, schema_editor):
    """ Move TranscriptVersion.data["cdot"] into its own column - once imports stopped rewriting unchanged rows
        (#2029) the key looked like the installed release when it's the one that row last changed in """
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")
    tv_table = TranscriptVersion._meta.db_table
    max_pk = TranscriptVersion.objects.aggregate(max_pk=Max("pk"))["max_pk"] or 0
    sql = f"""UPDATE {tv_table} SET modified_cdot_version = data->>'cdot', data = data - 'cdot'
              WHERE id >= %s AND id < %s AND data ? 'cdot'"""
    with schema_editor.connection.cursor() as cursor:
        for start in range(0, max_pk + 1, BATCH_SIZE):
            cursor.execute(sql, [start, start + BATCH_SIZE])


class Migration(migrations.Migration):
    # Commit each batch rather than holding row locks on ~2.2M TranscriptVersions in one transaction
    atomic = False

    dependencies = [
        ("genes", "0099_transcriptversion_modified_cdot_version"),
    ]

    operations = [
        migrations.RunPython(_move_data_cdot_to_modified_cdot_version, migrations.RunPython.noop),
    ]
