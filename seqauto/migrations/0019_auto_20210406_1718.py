# Generated by Django 3.1.3 on 2021-04-06 07:48

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0026_auto_20210406_1718'),
        ('snpdb', '0029_legacy_vcf_allele_frequency_percent'),
        ('seqauto', '0018_one_off_seqauto_variant_caller'),
    ]

    operations = [
        migrations.RenameField(
            model_name='goldcoveragesummary',
            old_name='original_transcript_id',
            new_name='original_transcript',
        ),
        migrations.AddField(
            model_name='goldcoveragesummary',
            name='transcript_version',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='genes.transcriptversion'),
        ),
        migrations.AlterField(
            model_name='enrichmentkit',
            name='canonical_transcript_collection',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='genes.canonicaltranscriptcollection'),
        ),
        migrations.AlterField(
            model_name='enrichmentkit',
            name='gene_list',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, to='genes.genelist'),
        ),
        migrations.AlterField(
            model_name='enrichmentkit',
            name='manufacturer',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='snpdb.manufacturer'),
        ),
    ]