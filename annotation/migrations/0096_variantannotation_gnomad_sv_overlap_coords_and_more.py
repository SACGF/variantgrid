# Generated by Django 4.2.10 on 2024-03-27 03:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0095_one_off_gnomad_sv_hash_check'),
    ]

    operations = [
        migrations.AddField(
            model_name='variantannotation',
            name='gnomad_sv_overlap_coords',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='annotationrun',
            name='pipeline_type',
            field=models.CharField(choices=[('S', 'Standard Short Variant'), ('C', 'Structural Variant')], default='S', max_length=1),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='pipeline_type',
            field=models.CharField(blank=True, choices=[('S', 'Standard Short Variant'), ('C', 'Structural Variant')], max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='vep_custom',
            field=models.CharField(choices=[('g', 'gnomAD2'), ('n', 'gnomAD3'), ('o', 'gnomAD4'), ('S', 'gnomAD_SV'), ('N', 'gnomAD_SV_name'), ('1', 'phastCons100way'), ('2', 'phastCons30way'), ('3', 'phastCons46way'), ('4', 'phyloP100way'), ('5', 'phyloP30way'), ('6', 'phyloP46way'), ('r', 'RepeatMasker'), ('t', 'TopMed'), ('u', 'UK10k'), ('c', 'COSMIC')], max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='columnvepfield',
            name='vep_plugin',
            field=models.CharField(choices=[('d', 'dbNSFP'), ('v', 'dbscSNV'), ('g', 'Grantham'), ('l', 'LoFtool'), ('n', 'Mastermind'), ('V', 'MaveDb'), ('m', 'MaxEntScan'), ('N', 'NMD'), ('a', 'SpliceAI'), ('s', 'SpliceRegion')], max_length=1, null=True),
        ),
    ]
