from django.db import migrations, models


def _flip_existing_vavs_false(apps, schema_editor):
    """Existing VAVs predate the optimisations, so default them to False until
    the matching `fix_historical_*` / `fix_columns_version*_damage_counts`
    management command flips them True per-partition."""
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.all().update(
        backfilled_spliceai_max_ds=False,
        backfilled_max_af=False,
        backfilled_damage_counts=False,
    )


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0150_drop_annotation_run_from_vgo_unique_together"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_spliceai_max_ds",
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_max_af",
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_damage_counts",
            field=models.BooleanField(default=True),
        ),
        migrations.RunPython(_flip_existing_vavs_false, migrations.RunPython.noop),
    ]
