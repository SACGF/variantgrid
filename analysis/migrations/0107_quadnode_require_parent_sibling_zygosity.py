from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0106_alter_populationnode_af_1kg_and_more"),
    ]

    operations = [
        migrations.RenameField(
            model_name="quadnode",
            old_name="require_zygosity",
            new_name="require_parent_zygosity",
        ),
        migrations.AddField(
            model_name="quadnode",
            name="require_sibling_zygosity",
            field=models.BooleanField(default=True),
        ),
    ]
