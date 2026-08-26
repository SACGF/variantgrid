from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0172_gene_fusion_variant_class"),
    ]

    operations = [
        migrations.AddField(
            model_name="annotationrun",
            name="dispatch_count",
            field=models.IntegerField(default=0),
        ),
    ]
