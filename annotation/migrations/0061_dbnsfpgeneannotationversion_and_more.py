# Generated by Django 4.0.6 on 2022-08-10 04:20

import django.db.models.deletion
import django_extensions.db.fields
import psqlextra.manager.manager
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0057_one_off_upgrade_pyhgvs'),
        ('annotation', '0060_one_off_gene_annotation_ontology_version'),
    ]

    operations = [
        migrations.CreateModel(
            name='DBNSFPGeneAnnotationVersion',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('version', models.TextField(primary_key=True, serialize=False)),
                ('md5_hash', models.CharField(max_length=32, unique=True)),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='variantannotation',
            name='loftool',
        ),
        migrations.RemoveField(
            model_name='varianttranscriptannotation',
            name='loftool',
        ),
        migrations.CreateModel(
            name='DBNSFPGeneAnnotation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('gene_damage_index_score', models.FloatField(null=True)),
                ('gene_damage_index_phred', models.FloatField(null=True)),
                ('phi', models.FloatField(null=True)),
                ('ghis', models.FloatField(null=True)),
                ('prec', models.FloatField(null=True)),
                ('hipred_score', models.FloatField(null=True)),
                ('gnomad_pli', models.FloatField(blank=True, null=True)),
                ('gnomad_prec', models.FloatField(blank=True, null=True)),
                ('gnomad_pnull', models.FloatField(blank=True, null=True)),
                ('loftool', models.FloatField(blank=True, null=True)),
                ('gene_indispensability_score', models.FloatField(blank=True, null=True)),
                ('hipred_prediction', models.BooleanField(null=True)),
                ('gene_indispensability_pred', models.BooleanField(null=True)),
                ('pathway_biocarta_full', models.TextField(blank=True, null=True)),
                ('pathway_consensus_pathdb', models.TextField(blank=True, null=True)),
                ('pathway_kegg_id', models.TextField(blank=True, null=True)),
                ('pathway_kegg_full', models.TextField(blank=True, null=True)),
                ('gwas_trait_association', models.TextField(blank=True, null=True)),
                ('go_biological_process', models.TextField(blank=True, null=True)),
                ('go_cellular_component', models.TextField(blank=True, null=True)),
                ('go_molecular_function', models.TextField(blank=True, null=True)),
                ('interactions_biogrid', models.TextField(blank=True, null=True)),
                ('interactions_consensus_pathdb', models.TextField(blank=True, null=True)),
                ('expression_egenetics', models.TextField(blank=True, null=True)),
                ('expression_gnf_atlas', models.TextField(blank=True, null=True)),
                ('essential_gene_crispr', models.CharField(blank=True, choices=[('E', 'Essential'), ('N', 'Non-essential phenotype-changing')], max_length=1, null=True)),
                ('essential_gene_crispr2', models.CharField(blank=True, choices=[('E', 'Essential'), ('N', 'Non-essential phenotype-changing'), ('S', 'Context-Specific essential')], max_length=1, null=True)),
                ('essential_gene_gene_trap', models.CharField(blank=True, choices=[('E', 'Essential'), ('N', 'Non-essential phenotype-changing'), ('H', 'HAP1-Specific essential'), ('K', 'KBM7-Specific essential')], max_length=1, null=True)),
                ('ensembl_transcript', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='ensembl_dbnsfp_gene', to='genes.transcript')),
                ('gene_symbol', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genes.genesymbol')),
                ('refseq_transcript', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='refseq_dbnsfp_gene', to='genes.transcript')),
                ('version', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='annotation.dbnsfpgeneannotationversion')),
            ],
            options={
                'abstract': False,
                'base_manager_name': 'objects',
            },
            managers=[
                ('objects', psqlextra.manager.manager.PostgresManager()),
            ],
        ),
        migrations.AddField(
            model_name='geneannotation',
            name='dbnsfp_gene',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='annotation.dbnsfpgeneannotation'),
        ),
        migrations.AddField(
            model_name='geneannotationversion',
            name='dbnsfp_gene_version',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='annotation.dbnsfpgeneannotationversion'),
        ),
    ]