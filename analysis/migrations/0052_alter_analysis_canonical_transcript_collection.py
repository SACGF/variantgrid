# Generated by Django 3.2.8 on 2021-12-06 04:05

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0048_one_off_upgrade_pyhgvs'),
        ('analysis', '0051_analysis_canonical_transcript_collection'),
    ]

    operations = [
        migrations.AlterField(
            model_name='analysis',
            name='canonical_transcript_collection',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.canonicaltranscriptcollection'),
        ),
    ]