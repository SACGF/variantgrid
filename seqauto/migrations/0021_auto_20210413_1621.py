# Generated by Django 3.1.3 on 2021-04-13 06:51

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('seqauto', '0020_remove_enrichmentkit_bed_file'),
    ]

    operations = [
        migrations.RenameField(
            model_name='qcexecsummary',
            old_name='percent_10x',
            new_name='percent_10x_goi',
        ),
        migrations.RenameField(
            model_name='qcexecsummary',
            old_name='percent_20x',
            new_name='percent_20x_goi',
        ),
        migrations.RenameField(
            model_name='qcexecsummary',
            old_name='percent_250x',
            new_name='percent_250x_goi',
        ),
        migrations.RenameField(
            model_name='qcexecsummary',
            old_name='percent_500x',
            new_name='percent_500x_goi',
        ),
        migrations.RenameField(
            model_name='qcexecsummary',
            old_name='duplicated_alignable_reads',
            new_name='percent_duplication',
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='deduplicated_reads',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='percent_error_rate',
            field=models.FloatField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='percent_map_to_diff_chr',
            field=models.FloatField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='percent_reads',
            field=models.FloatField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='percent_softclip',
            field=models.FloatField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='reads',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='sample_id_lod',
            field=models.FloatField(null=True),
        ),
        migrations.AddField(
            model_name='qcexecsummary',
            name='sex_match',
            field=models.TextField(null=True),
        ),
    ]
