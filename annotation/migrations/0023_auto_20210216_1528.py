# Generated by Django 3.1 on 2021-02-16 04:58

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0019_new_columns_gnomad3'),
        ('annotation', '0022_one_off_manual_tasks_gene_annotation'),
    ]

    operations = [
        migrations.AddField(
            model_name='columnvepfield',
            name='genome_build',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='snpdb.genomebuild'),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad2_liftover_af',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_ac',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_an',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_popmax_ac',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_popmax_an',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_popmax_hom_alt',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='variant_grid_column',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='snpdb.variantgridcolumn'),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='vep_custom',
            field=models.CharField(choices=[('g', 'gnomAD2'), ('n', 'gnomAD3'), ('1', 'phastCons100way'), ('2', 'phastCons30way'), ('3', 'phastCons46way'), ('4', 'phyloP100way'), ('5', 'phyloP30way'), ('6', 'phyloP46way'), ('r', 'RepeatMasker'), ('t', 'TopMed'), ('u', 'UK10k'), ('c', 'COSMIC')], max_length=1, null=True),
        ),
    ]
