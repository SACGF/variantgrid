# Generated by Django 4.2.15 on 2024-12-05 01:42

from django.db import migrations


def _vep_column_sift_37_38_only(apps, _schema_editor):
    ColumnVEPField = apps.get_model("annotation", "ColumnVEPField")
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")

    GENOME_BUILDS_37_38_ONLY = [
        'Aloft_Confidence',
        'Aloft_pred',
        'Aloft_prob_Dominant',
        'Aloft_prob_Recessive',
        'Aloft_prob_Tolerant',
        'AlphaMissense_rankscore',
        'BayesDel_noAF_rankscore',
        'CADD_raw_rankscore',
        'ClinPred_rankscore',
        'Ensembl_transcriptid',
        'GERP++_RS',
        'Interpro_domain',
        'Mastermind_MMID3',
        'Mastermind_counts',
        'MetaLR_rankscore',
        'REVEL_rankscore',
        'SIFT',
        'SpliceAI_pred_DP_AG',
        'SpliceAI_pred_DP_AL',
        'SpliceAI_pred_DP_DG',
        'SpliceAI_pred_DP_DL',
        'SpliceAI_pred_DS_AG',
        'SpliceAI_pred_DS_AL',
        'SpliceAI_pred_DS_DG',
        'SpliceAI_pred_DS_DL',
        'SpliceAI_pred_SYMBOL',
        'VEST4_rankscore',
        'ada_score',
        'rf_score',
    ]

    genome_build_37 = GenomeBuild.objects.get(pk="GRCh37")
    genome_build_38 = GenomeBuild.objects.get(pk="GRCh38")
    for source_field in GENOME_BUILDS_37_38_ONLY:
        for cvf in ColumnVEPField.objects.filter(source_field=source_field):
            original_column = cvf.column

            # Make original one 37
            cvf.column = original_column + "_37"
            cvf.genome_build = genome_build_37
            cvf.save()

            cvf.column = original_column + "_38"
            cvf.pk = None  # Make new 38
            cvf.genome_build = genome_build_38
            cvf.save()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0120_alter_variantannotationversion_thousand_genomes'),
    ]

    operations = [
        migrations.RunPython(_vep_column_sift_37_38_only)
    ]