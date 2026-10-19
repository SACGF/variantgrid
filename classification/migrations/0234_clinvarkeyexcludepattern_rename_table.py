from django.db import migrations


class Migration(migrations.Migration):
    """ #2276 - rename snpdb_clinvarkeyexcludepattern to classification_clinvarkeyexcludepattern. A rename only: its
        constraint/index names and id sequence keep their snpdb_ prefix """

    dependencies = [
        ("classification", "0233_clinvarkeyexcludepattern"),
    ]

    operations = [
        migrations.AlterModelTable(
            name="clinvarkeyexcludepattern",
            table=None,
        ),
        migrations.AlterModelOptions(
            name="clinvarkeyexcludepattern",
            options={"get_latest_by": "modified"},
        ),
    ]
