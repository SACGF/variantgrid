# Generated by Django 4.0.4 on 2022-07-06 03:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0077_rename_continued_discordance_text_discordancereport_notes'),
    ]

    operations = [
        migrations.AddField(
            model_name='classificationimportrun',
            name='missing_row_count',
            field=models.IntegerField(blank=True, null=True),
        ),
    ]