# Generated by Django 3.2.8 on 2021-12-01 06:42

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0048_one_off_upgrade_pyhgvs'),
        ('analysis', '0050_auto_20210622_1459'),
    ]

    operations = [
        migrations.AddField(
            model_name='analysis',
            name='canonical_transcript_collection',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.canonicaltranscriptcollection'),
        ),
    ]
