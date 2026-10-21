from django.db import migrations


class Migration(migrations.Migration):
    """ #2278 - rename annotation_cachedwebresource to snpdb_cachedwebresource. A rename only: the FKs
        referencing it follow the table, and its constraint/index names keep their annotation_ prefix """

    dependencies = [
        ("snpdb", "0289_cachedwebresource"),
        ("annotation", "0198_cached_web_resource_to_snpdb"),
    ]

    operations = [
        migrations.AlterModelTable(
            name="cachedwebresource",
            table=None,
        ),
    ]
