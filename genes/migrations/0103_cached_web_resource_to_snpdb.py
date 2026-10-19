import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    """ #2278 - repoint FKs at snpdb.CachedWebResource. State only: the table and constraints are unchanged """

    dependencies = [
        ("genes", "0102_spliceevent_civic_variant_id"),
        ("snpdb", "0289_cachedwebresource"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.AlterField(
                    model_name=model_name,
                    name="cached_web_resource",
                    field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="snpdb.cachedwebresource"),
                )
                for model_name in ("cachedthirdpartygenelist", "gnomadgeneconstraint", "lrgrefseqgene", "uniprot")
            ],
        ),
    ]
