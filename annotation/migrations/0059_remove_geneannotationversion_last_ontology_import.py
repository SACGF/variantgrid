# Generated by Django 4.0.6 on 2022-07-28 04:48

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0058_one_off_gene_annotation_assign_ontology_version'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='geneannotationversion',
            name='last_ontology_import',
        ),
    ]
