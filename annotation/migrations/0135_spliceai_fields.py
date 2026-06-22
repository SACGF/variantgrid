from django.db import migrations, models


def _set_existing_vavs_to_raw(apps, schema_editor):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.filter(spliceai__isnull=True).update(spliceai="raw 1.3")


def _noop(apps, schema_editor):
    pass


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0134_annotsv_fields"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="spliceai",
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="variantannotation",
            name="spliceai_max_ds",
            field=models.FloatField(blank=True, db_index=True, null=True),
        ),
        migrations.RunPython(_set_existing_vavs_to_raw, _noop),
    ]
