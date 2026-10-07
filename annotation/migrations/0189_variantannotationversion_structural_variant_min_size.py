from django.db import migrations, models


class Migration(migrations.Migration):
    """ Existing versions annotated every symbolic variant as structural, and every symbolic variant was
        1000bp or more (#1358) """

    dependencies = [
        ("annotation", "0188_one_off_fix_annotation_vep_too_long"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="structural_variant_min_size",
            field=models.IntegerField(default=1000),
            preserve_default=False,
        ),
    ]
