# Generated by Django 3.1.6 on 2021-07-30 03:22

import django.utils.timezone
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0043_auto_20210723_1441'),
    ]

    operations = [
        migrations.AddField(
            model_name='clinvarexport',
            name='last_evaluated',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
        migrations.AddField(
            model_name='clinvarexport',
            name='submission_body_validated',
            field=models.JSONField(default=dict),
        ),
    ]
