# Generated by Django 4.0.6 on 2022-07-27 05:44

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0078_alter_vcf_options'),
    ]

    operations = [
        migrations.RenameField(
            model_name='allelemergelog',
            old_name='conversion_tool',
            new_name='allele_linking_tool',
        ),
        migrations.RenameField(
            model_name='variantallele',
            old_name='conversion_tool',
            new_name='allele_linking_tool',
        ),
    ]
