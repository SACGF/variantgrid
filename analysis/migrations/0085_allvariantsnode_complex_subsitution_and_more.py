# Generated by Django 4.2.10 on 2024-05-10 07:54

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0084_remove_karyomappinggene_gene_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='allvariantsnode',
            name='complex_subsitution',
            field=models.BooleanField(blank=True, default=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='indels',
            field=models.BooleanField(blank=True, default=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='snps',
            field=models.BooleanField(blank=True, default=True),
        ),
        migrations.AddField(
            model_name='allvariantsnode',
            name='structural_variants',
            field=models.BooleanField(blank=True, default=True),
        ),
    ]