# Generated by Django 4.2.2 on 2023-11-24 04:41

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0079_variantannotation_alphamissense_class_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='columnvepfield',
            name='pipeline_type',
            field=models.CharField(choices=[('S', 'Standard Short Variant'), ('C', 'CNV')], default='S', max_length=1),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='faf95',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='faf99',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='fafmax_faf95_max',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotation',
            name='fafmax_faf99_max',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='vep_custom',
            field=models.CharField(choices=[('g', 'gnomAD2'), ('n', 'gnomAD3'), ('o', 'gnomAD4'), ('1', 'phastCons100way'), ('2', 'phastCons30way'), ('3', 'phastCons46way'), ('4', 'phyloP100way'), ('5', 'phyloP30way'), ('6', 'phyloP46way'), ('r', 'RepeatMasker'), ('t', 'TopMed'), ('u', 'UK10k'), ('c', 'COSMIC')], max_length=1, null=True),
        ),
    ]