from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0178_one_off_recalculate_symbolic_hgvs"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotation",
            name="open_targets_gwas_l2g_scores",
            field=models.TextField(blank=True, null=True),
        ),
    ]
