# Generated by Django 4.0.2 on 2022-02-23 02:00

import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0058_importedwikicollection_importedwiki'),
    ]

    operations = [
        migrations.CreateModel(
            name='CohortGenotypeCommonFilterVersion',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('gnomad_version', models.TextField()),
                ('gnomad_af_min', models.FloatField()),
                ('clinical_significance_max', models.CharField(max_length=1, null=True)),
                ('genome_build', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.genomebuild')),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='cohortgenotypecollection',
            name='common_collection',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='uncommon', to='snpdb.cohortgenotypecollection'),
        ),
        migrations.AddField(
            model_name='cohortgenotypecollection',
            name='common_filter',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.PROTECT, to='snpdb.cohortgenotypecommonfilterversion'),
        ),
        migrations.CreateModel(
            name='CommonVariantClassified',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('common_filter', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.cohortgenotypecommonfilterversion')),
                ('variant', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.variant')),
            ],
            options={
                'unique_together': {('variant', 'common_filter')},
            },
        ),
    ]