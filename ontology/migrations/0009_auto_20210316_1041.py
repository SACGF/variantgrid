# Generated by Django 3.1.6 on 2021-03-16 00:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ontology', '0008_load_ontology'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ontologyterm',
            name='ontology_service',
            field=models.CharField(choices=[('MONDO', 'MONDO'), ('OMIM', 'OMIM'), ('HP', 'HP'), ('HGNC', 'HGNC'), ('DOID', 'DOID'), ('Orphanet', 'Orphanet')], max_length=10),
        ),
    ]