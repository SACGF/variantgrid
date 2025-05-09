# Generated by Django 4.2.18 on 2025-02-03 04:43

import django.contrib.postgres.fields
import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0161_alter_cohortgenotype_unique_together'),
    ]

    operations = [
        migrations.CreateModel(
            name='VCFLengthStatsCollection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('code_version', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.samplestatscodeversion')),
                ('vcf', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='snpdb.vcf')),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='VCFLengthStats',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('variant_class', models.CharField(choices=[('SN', 'SNV'), ('GM', 'genetic_marker'), ('SU', 'substitution'), ('TR', 'tandem_repeat'), ('AI', 'Alu_insertion'), ('CA', 'complex_structural_alteration'), ('CS', 'complex_substitution'), ('CG', 'copy_number_gain'), ('CL', 'copy_number_loss'), ('CN', 'copy_number_variation'), ('DU', 'duplication'), ('IB', 'interchromosomal_breakpoint'), ('IT', 'interchromosomal_translocation'), ('CB', 'intrachromosomal_breakpoint'), ('CT', 'intrachromosomal_translocation'), ('IN', 'inversion'), ('LO', 'loss_of_heterozygosity'), ('MD', 'mobile_element_deletion'), ('MI', 'mobile_element_insertion'), ('NI', 'novel_sequence_insertion'), ('ST', 'short_tandem_repeat_variation'), ('TD', 'tandem_duplication'), ('TL', 'translocation'), ('DE', 'deletion'), ('ND', 'indel'), ('IS', 'insertion'), ('SA', 'sequence_alteration'), ('PR', 'probe')], max_length=2, null=True)),
                ('is_log', models.BooleanField(default=False)),
                ('histogram_counts', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(), size=None)),
                ('histogram_bin_edges', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(), size=None)),
                ('collection', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.vcflengthstatscollection')),
            ],
            options={
                'unique_together': {('collection', 'variant_class')},
            },
        ),
    ]
