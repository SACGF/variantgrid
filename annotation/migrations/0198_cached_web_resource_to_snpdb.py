import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    """ #2278 - CachedWebResource moved to snpdb. State only: snpdb's model keeps the annotation_cachedwebresource
        table, so FKs are repointed and the model deleted without touching the database """

    dependencies = [
        ("annotation", "0197_variantannotationversion_vep_version_mismatch"),
        ("genes", "0103_cached_web_resource_to_snpdb"),
        ("snpdb", "0289_cachedwebresource"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.AlterField(
                    model_name="clinvarcitationscollection",
                    name="cached_web_resource",
                    field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to="snpdb.cachedwebresource"),
                ),
                migrations.AlterField(
                    model_name="genepubmedcount",
                    name="cached_web_resource",
                    field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="snpdb.cachedwebresource"),
                ),
                migrations.DeleteModel(
                    name="CachedWebResource",
                ),
            ],
        ),
    ]
