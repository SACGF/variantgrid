# Generated by Django 4.0.3 on 2022-04-13 02:55

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0053_alter_canonicaltranscript_transcript_and_more'),
        ('annotation', '0042_one_off_variant_annotation_version_gnomad'),
    ]

    operations = [
        migrations.AlterField(
            model_name='variantannotation',
            name='transcript',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.transcript'),
        ),
        migrations.AlterField(
            model_name='variantannotation',
            name='transcript_version',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.transcriptversion'),
        ),
        migrations.AlterField(
            model_name='varianttranscriptannotation',
            name='transcript',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.transcript'),
        ),
        migrations.AlterField(
            model_name='varianttranscriptannotation',
            name='transcript_version',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.transcriptversion'),
        ),
    ]
