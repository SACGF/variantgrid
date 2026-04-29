from django.db import migrations, models


_TD = [('T', 'Tolerated'), ('D', 'Damaging')]
_ALPHAMISSENSE = [('b', 'likely_benign'), ('a', 'ambiguous'), ('p', 'likely_pathogenic')]


def _new_fields():
    return [
        ('alphamissense_pred', models.CharField(blank=True, choices=_ALPHAMISSENSE, max_length=1, null=True)),
        ('alphamissense_score', models.FloatField(blank=True, null=True)),
        ('bayesdel_noaf_score', models.FloatField(blank=True, null=True)),
        ('cadd_raw', models.FloatField(blank=True, null=True)),
        ('clinpred_pred', models.CharField(blank=True, choices=_TD, max_length=1, null=True)),
        ('clinpred_score', models.FloatField(blank=True, null=True)),
        ('metarnn_pred', models.CharField(blank=True, choices=_TD, max_length=1, null=True)),
        ('metarnn_score', models.FloatField(blank=True, null=True)),
        ('mpc_score', models.FloatField(blank=True, null=True)),
        ('mutpred2_score', models.FloatField(blank=True, null=True)),
        ('mutpred2_top5_mechanisms', models.TextField(blank=True, null=True)),
        ('primateai_pred', models.CharField(blank=True, choices=_TD, max_length=1, null=True)),
        ('primateai_score', models.FloatField(blank=True, null=True)),
        ('varity_er_score', models.FloatField(blank=True, null=True)),
        ('varity_r_score', models.FloatField(blank=True, null=True)),
        ('vest4_score', models.FloatField(blank=True, null=True)),
    ]


def _add_field_ops():
    ops = []
    for model_name in ('variantannotation', 'varianttranscriptannotation'):
        for name, field in _new_fields():
            ops.append(migrations.AddField(model_name=model_name, name=name, field=field))
    return ops


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0130_one_off_fix_annotation_sv_overlaps'),
    ]

    operations = _add_field_ops()
