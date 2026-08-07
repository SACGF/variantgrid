from django.db import migrations


def _flip_columns_version3_true(apps, schema_editor):
    """columns_version=3 VAVs never needed a damage-counts backfill: their
    predictions_num_pathogenic / predictions_num_benign aggregates were populated
    at insert time (the six shared rankscores plus an alphamissense contribution),
    and unlike v2 (0054) and v4 (0146) there was never a one-off recompute migration
    for v3. Migration 0151 flipped every VAV's backfilled_damage_counts False as a
    blanket conservative default, which incorrectly disabled the DamageNode
    damage_predictions_min filter for v3 and pointed at a `fix_columns_version3_*`
    command that was never written. Their data is already correct, so just re-enable
    the flag - no recompute required."""
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.filter(columns_version=3).update(backfilled_damage_counts=True)


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0154_annotationrun_attempt_count_and_more"),
    ]

    operations = [
        migrations.RunPython(_flip_columns_version3_true, migrations.RunPython.noop),
    ]
