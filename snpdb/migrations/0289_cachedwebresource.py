import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):
    """ #2278 - CachedWebResource moves from annotation to snpdb in state only, keeping its table """

    dependencies = [
        ("snpdb", "0288_one_off_allele_liftover_stuck_processing"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.CreateModel(
                    name="CachedWebResource",
                    fields=[
                        ("created", django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name="created")),
                        ("modified", django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name="modified")),
                        ("name", models.TextField(primary_key=True, serialize=False)),
                        ("description", models.TextField(blank=True)),
                        ("import_status", models.CharField(choices=[("C", "created"), ("I", "importing"), ("R", "Requires user input"), ("E", "error"), ("S", "success"), ("M", "Marked For Deletion"), ("D", "Deleting")], default="C", max_length=1)),
                    ],
                    options={
                        "db_table": "annotation_cachedwebresource",
                        "get_latest_by": "modified",
                        "abstract": False,
                    },
                ),
            ],
        ),
    ]
