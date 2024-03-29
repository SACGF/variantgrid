# Generated by Django 3.2.6 on 2021-09-14 12:08

import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0038_lrgrefseqgene'),
    ]

    operations = [
        migrations.CreateModel(
            name='TranscriptVersionSequenceInfoFastaFileImport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('md5_hash', models.CharField(max_length=32, unique=True)),
                ('annotation_consortium', models.CharField(choices=[('R', 'RefSeq'), ('E', 'Ensembl')], max_length=1)),
                ('filename', models.TextField()),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='TranscriptVersionSequenceInfo',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('version', models.IntegerField()),
                ('api_response', models.TextField(null=True)),
                ('sequence', models.TextField()),
                ('length', models.IntegerField()),
                ('fasta_import', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='genes.transcriptversionsequenceinfofastafileimport')),
                ('transcript', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genes.transcript')),
            ],
            options={
                'unique_together': {('transcript', 'version')},
            },
        ),
    ]
