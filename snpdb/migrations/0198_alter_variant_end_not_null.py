from django.db import migrations, models


def _require_variant_end_populated(apps, _schema_editor):
    Variant = apps.get_model("snpdb", "Variant")
    if Variant.objects.filter(end__isnull=True).exists():
        raise RuntimeError(
            "Cannot make Variant.end NOT NULL - some variants still have end=NULL.\n"
            "Populate it first (this can take a long time, so run it outside migrations in screen):\n"
            "    python3 manage.py one_off_fix_variant_end\n"
            "then re-run migrate."
        )


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0197_allvariantsfilter'),
    ]

    operations = [
        migrations.RunPython(_require_variant_end_populated, migrations.RunPython.noop),
        migrations.AlterField(
            model_name='variant',
            name='end',
            field=models.IntegerField(),
        ),
    ]
