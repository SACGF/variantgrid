# Generated by Django 3.1 on 2021-01-12 04:55

import django.utils.timezone
import model_utils.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('genes', '0009_one_off_fix_unknown_genes'),
    ]

    operations = [
        migrations.CreateModel(
            name='OntologyImport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('processed_date', models.DateTimeField(auto_created=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                ('ontology_service', models.CharField(choices=[('MONDO', 'MONDO'), ('OMIM', 'OMIM'), ('HPO', 'HPO'), ('PAAU', 'PanelApp AU')], max_length=5)),
                ('filename', models.TextField()),
                ('context', models.TextField()),
                ('hash', models.TextField()),
                ('completed', models.BooleanField(default=False)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='OntologyTerm',
            fields=[
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                ('id', models.TextField(primary_key=True, serialize=False)),
                ('ontology_service', models.CharField(choices=[('MONDO', 'MONDO'), ('OMIM', 'OMIM'), ('HPO', 'HPO'), ('PAAU', 'PanelApp AU')], max_length=5)),
                ('index', models.IntegerField()),
                ('name', models.TextField(blank=True, null=True)),
                ('definition', models.TextField(blank=True, null=True)),
                ('extra', models.JSONField(blank=True, null=True)),
                ('from_import', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='ontology.ontologyimport')),
            ],
            options={
                'unique_together': {('ontology_service', 'index')},
            },
        ),
        migrations.CreateModel(
            name='OntologyTermRelation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                ('relation', models.CharField(choices=[('is_a', 'is a'), ('exact', 'exact'), ('close', 'close'), ('broad', 'broad'), ('narrow', 'narrow')], max_length=10)),
                ('dest_term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ontology.ontologyterm')),
                ('from_import', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='ontology.ontologyimport')),
                ('source_term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='subject', to='ontology.ontologyterm')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='OntologyTermGeneRelation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                ('relation', models.TextField()),
                ('extra', models.JSONField(blank=True, null=True)),
                ('deletion_date', models.DateTimeField(blank=True, null=True)),
                ('from_import', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='ontology.ontologyimport')),
                ('gene_symbol', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genes.genesymbol')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ontology.ontologyterm')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
