from django.db import migrations


class Migration(migrations.Migration):
    """ #2276 - rename snpdb_clinvarkey / snpdb_clinvarkeyexcludepattern to classification_. A rename only: the FKs
        referencing them follow the tables, and their constraint/index names keep the snpdb_ prefix """

    dependencies = [
        ("classification", "0233_classification_settings"),
        ("snpdb", "0291_move_classification_settings"),
    ]

    operations = [
        migrations.AlterModelTable(
            name="clinvarkey",
            table=None,
        ),
        migrations.AlterModelTable(
            name="clinvarkeyexcludepattern",
            table=None,
        ),
        migrations.AlterModelOptions(
            name="clinvarkeyexcludepattern",
            options={"get_latest_by": "modified"},
        ),
    ]
