# Generated by Django 4.2.5 on 2023-11-24 06:08

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0120_once_off_some_somatic_keys'),
    ]

    operations = [
        migrations.AddField(
            model_name='clinvarallele',
            name='allele_origin_bucket',
            field=models.CharField(choices=[('G', 'Germline'), ('S', 'Somatic'), ('U', 'Unknown')], default='G', max_length=1),
        ),
        migrations.AddField(
            model_name='clinvarexport',
            name='allele_origin_bucket',
            field=models.CharField(choices=[('G', 'Germline'), ('S', 'Somatic'), ('U', 'Unknown')], default='G', max_length=1),
        ),
    ]