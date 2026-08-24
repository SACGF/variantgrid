from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0214_vcfsourcesettings_genome_build"),
    ]

    operations = [
        migrations.RenameField(
            model_name="jobscontrol",
            old_name="last_autopause_boot_time",
            new_name="last_boot_time",
        ),
        migrations.AddField(
            model_name="jobscontrol",
            name="previous_boot_time",
            field=models.BigIntegerField(blank=True, null=True),
        ),
    ]
