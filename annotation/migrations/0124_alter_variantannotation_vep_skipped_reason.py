# Generated by Django 4.2.18 on 2025-01-28 23:32

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0123_one_off_fix_annotation_sv_c_hgvs'),
    ]

    operations = [
        migrations.AlterField(
            model_name='variantannotation',
            name='vep_skipped_reason',
            field=models.CharField(blank=True, choices=[('c', 'Unknown Contig'), ('i', 'Incomplete'), ('u', 'Unknown'), ('l', 'Too Long')], max_length=1, null=True),
        ),
    ]
