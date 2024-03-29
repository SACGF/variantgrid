# Generated by Django 3.2.1 on 2021-05-24 02:55

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('seqauto', '0025_auto_20210429_1030'),
    ]

    operations = [
        migrations.AddField(
            model_name='enrichmentkit',
            name='sample_variants_type',
            field=models.CharField(choices=[('U', 'Unknown'), ('G', 'Germline'), ('M', 'Mixed (Single Sample)'), ('S', 'Somatic only (Tumor minus normal)')], default='U', max_length=1),
        ),
        migrations.AddField(
            model_name='enrichmentkit',
            name='variant_zygosity_count',
            field=models.BooleanField(default=True),
        ),
    ]
