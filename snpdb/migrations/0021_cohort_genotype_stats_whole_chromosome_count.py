from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0020_sample_ploidy"),
    ]

    operations = [
        migrations.AddField(
            model_name="cohortgenotypestats",
            name="whole_chromosome_count",
            field=models.IntegerField(default=0),
        ),
    ]
