# Generated by Django 4.0.7 on 2023-05-17 05:24

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('review', '0001_initial'),
        ('classification', '0102_one_off_update_ekey_descriptions'),
    ]

    operations = [
        migrations.AddField(
            model_name='discordancereport',
            name='reviews',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='review.reviewedobject'),
        ),
    ]
