# Generated by Django 3.1 on 2020-11-18 04:40

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0003_transcript_version_json_del_gene_name'),
        ('analysis', '0006_one_off_fix_analysis_permissions'),
    ]

    operations = [
        migrations.AddField(
            model_name='allvariantsnode',
            name='gene_symbol',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='genes.genesymbol'),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='max_het_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='max_hom_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='max_ref_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='maximum_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='min_het_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='min_hom_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='min_ref_count',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='minimum_count',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='reference',
            field=models.BooleanField(blank=True, default=False),
        ),
        migrations.AlterField(
            model_name='cohortnode',
            name='maximum_count',
            field=models.IntegerField(blank=True, null=True),
        ),
    ]
