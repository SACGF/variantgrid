from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0179_open_targets_gwas_l2g_scores"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotation",
            name="open_targets_is_lead",
            field=models.TextField(blank=True, null=True),
        ),
    ]
