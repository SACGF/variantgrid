from django.db import migrations, models


def _new_variant_annotation_fields():
    return [
        # ProtVar
        ('protvar_stability', models.FloatField(blank=True, null=True)),
        ('protvar_pocket', models.TextField(blank=True, null=True)),
        ('protvar_int', models.TextField(blank=True, null=True)),
        # Open Targets
        ('open_targets_gwas_l2g_score', models.FloatField(blank=True, null=True)),
        ('open_targets_gwas_gene_id', models.TextField(blank=True, null=True)),
        ('open_targets_gwas_diseases', models.TextField(blank=True, null=True)),
        ('open_targets_study_type', models.TextField(blank=True, null=True)),
        ('open_targets_study_id', models.TextField(blank=True, null=True)),
        ('open_targets_variant_id', models.TextField(blank=True, null=True)),
        ('open_targets_qtl_gene_id', models.TextField(blank=True, null=True)),
        ('open_targets_qtl_biosample', models.TextField(blank=True, null=True)),
        # EVE / popEVE
        ('eve_score', models.FloatField(blank=True, null=True)),
        ('eve_class', models.TextField(blank=True, null=True)),
        ('popeve_score', models.FloatField(blank=True, null=True)),
        # PromoterAI
        ('promoter_ai_score', models.FloatField(blank=True, null=True)),
        ('promoter_ai_tss_pos', models.IntegerField(blank=True, null=True)),
    ]


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0155_fix_columns_version3_backfilled_damage_counts'),
    ]

    operations = [
        migrations.AddField(
            model_name='variantannotationversion',
            name='open_targets',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantannotationversion',
            name='popeve',
            field=models.TextField(blank=True, null=True),
        ),
        *[
            migrations.AddField(model_name='variantannotation', name=name, field=field)
            for name, field in _new_variant_annotation_fields()
        ],
    ]
