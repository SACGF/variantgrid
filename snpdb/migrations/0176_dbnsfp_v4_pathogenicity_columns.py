from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data


_TRANSCRIPT = 'T'

_NEW_COLUMNS = [
    {'grid_column_name': 'alphamissense_pred',
     'variant_column': 'variantannotation__alphamissense_pred',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'AlphaMissense pred',
     'description': '<a href="https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#alphamissense">AlphaMissense</a> categorical prediction (likely_benign / ambiguous / likely_pathogenic). From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'alphamissense_score',
     'variant_column': 'variantannotation__alphamissense_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'AlphaMissense score',
     'description': 'AlphaMissense raw pathogenicity score, range 0 to 1. Higher = more likely pathogenic. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'bayesdel_noaf_score',
     'variant_column': 'variantannotation__bayesdel_noaf_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'BayesDel noAF score',
     'description': 'BayesDel meta-score (no allele frequency variant), used for predicting deleteriousness of missense variants. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'cadd_raw',
     'variant_column': 'variantannotation__cadd_raw',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'CADD raw',
     'description': 'Raw CADD C-score (uncalibrated). See Kircher et al. (2014) Nature Genetics 46(3):310-5. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'clinpred_pred',
     'variant_column': 'variantannotation__clinpred_pred',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'ClinPred pred',
     'description': 'ClinPred categorical prediction: D=Damaging, T=Tolerated. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'clinpred_score',
     'variant_column': 'variantannotation__clinpred_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'ClinPred score',
     'description': 'ClinPred raw score, range 0 to 1. Higher = more likely pathogenic. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'metarnn_pred',
     'variant_column': 'variantannotation__metarnn_pred',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'MetaRNN pred',
     'description': 'MetaRNN categorical prediction: D=Damaging, T=Tolerated. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'metarnn_score',
     'variant_column': 'variantannotation__metarnn_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'MetaRNN score',
     'description': 'MetaRNN raw score, range 0 to 1. Higher = more likely pathogenic. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'mpc_score',
     'variant_column': 'variantannotation__mpc_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'MPC score',
     'description': 'MPC (Missense badness, PolyPhen-2, and Constraint) score. Higher = more deleterious. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'mutpred2_score',
     'variant_column': 'variantannotation__mutpred2_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'MutPred2 score',
     'description': 'MutPred2 general pathogenicity score, range 0 to 1. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'mutpred2_top5_mechanisms',
     'variant_column': 'variantannotation__mutpred2_top5_mechanisms',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'MutPred2 top5 mechanisms',
     'description': 'MutPred2 top 5 predicted molecular mechanisms of pathogenicity (semicolon-separated). From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'primateai_pred',
     'variant_column': 'variantannotation__primateai_pred',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'PrimateAI pred',
     'description': 'PrimateAI categorical prediction: D=Damaging, T=Tolerated. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'primateai_score',
     'variant_column': 'variantannotation__primateai_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'PrimateAI score',
     'description': 'PrimateAI raw score, range 0 to 1. Higher = more likely pathogenic. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'varity_er_score',
     'variant_column': 'variantannotation__varity_er_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'VARITY ER score',
     'description': 'VARITY_ER score (rare exomes training set) for missense pathogenicity, range 0 to 1. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'varity_r_score',
     'variant_column': 'variantannotation__varity_r_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'VARITY R score',
     'description': 'VARITY_R score (rare training set) for missense pathogenicity, range 0 to 1. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'vest4_score',
     'variant_column': 'variantannotation__vest4_score',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'VEST4 score',
     'description': 'VEST4 score for missense pathogenicity, range 0 to 1. Higher = more likely pathogenic. From dbNSFP.',
     'model_field': True,
     'queryset_field': True},
]

_NEW_COLUMN_VCF_INFO = [
    {'info_id': 'AlphaMissense_pred',     'column_id': 'alphamissense_pred',       'number': 1,    'type': 'S', 'description': 'AlphaMissense categorical prediction'},
    {'info_id': 'AlphaMissense_score',    'column_id': 'alphamissense_score',      'number': 1,    'type': 'F', 'description': 'AlphaMissense raw pathogenicity score'},
    {'info_id': 'BayesDel_noAF_score',    'column_id': 'bayesdel_noaf_score',      'number': 1,    'type': 'F', 'description': 'BayesDel score (no allele frequency)'},
    {'info_id': 'CADD_raw',               'column_id': 'cadd_raw',                 'number': 1,    'type': 'F', 'description': 'CADD raw (uncalibrated) C-score'},
    {'info_id': 'ClinPred_pred',          'column_id': 'clinpred_pred',            'number': 1,    'type': 'S', 'description': 'ClinPred categorical prediction (D/T)'},
    {'info_id': 'ClinPred_score',         'column_id': 'clinpred_score',           'number': 1,    'type': 'F', 'description': 'ClinPred raw score'},
    {'info_id': 'MetaRNN_pred',           'column_id': 'metarnn_pred',             'number': 1,    'type': 'S', 'description': 'MetaRNN categorical prediction (D/T)'},
    {'info_id': 'MetaRNN_score',          'column_id': 'metarnn_score',            'number': 1,    'type': 'F', 'description': 'MetaRNN raw score'},
    {'info_id': 'MPC_score',              'column_id': 'mpc_score',                'number': 1,    'type': 'F', 'description': 'MPC (Missense badness, PolyPhen-2, Constraint) score'},
    {'info_id': 'MutPred2_score',         'column_id': 'mutpred2_score',           'number': 1,    'type': 'F', 'description': 'MutPred2 general score'},
    {'info_id': 'MutPred2_top5_mechanisms', 'column_id': 'mutpred2_top5_mechanisms', 'number': None, 'type': 'S', 'description': 'MutPred2 top 5 predicted molecular mechanisms'},
    {'info_id': 'PrimateAI_pred',         'column_id': 'primateai_pred',           'number': 1,    'type': 'S', 'description': 'PrimateAI categorical prediction (D/T)'},
    {'info_id': 'PrimateAI_score',        'column_id': 'primateai_score',          'number': 1,    'type': 'F', 'description': 'PrimateAI raw score'},
    {'info_id': 'VARITY_ER_score',        'column_id': 'varity_er_score',          'number': 1,    'type': 'F', 'description': 'VARITY_ER score (rare exomes set)'},
    {'info_id': 'VARITY_R_score',         'column_id': 'varity_r_score',           'number': 1,    'type': 'F', 'description': 'VARITY_R score (rare set)'},
    {'info_id': 'VEST4_score',            'column_id': 'vest4_score',              'number': 1,    'type': 'F', 'description': 'VEST4 score'},
]


def _new_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COLUMNS)])
    bulk_insert_class_data(apps, "snpdb", [("ColumnVCFInfo", _NEW_COLUMN_VCF_INFO)])


def _reverse_new_columns(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    info_ids = [c['info_id'] for c in _NEW_COLUMN_VCF_INFO]
    apps.get_model("snpdb", "ColumnVCFInfo").objects.filter(info_id__in=info_ids).delete()
    apps.get_model("snpdb", "VariantGridColumn").objects.filter(grid_column_name__in=grid_names).delete()


def _append_to_all_columns(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    all_columns = CustomColumnsCollection.objects.get(name='All columns')
    sort_order_max = all_columns.customcolumn_set.aggregate(Max("sort_order"))["sort_order__max"] or 0
    new = []
    for i, c in enumerate(_NEW_COLUMNS):
        new.append(CustomColumn(custom_columns_collection=all_columns,
                                sort_order=sort_order_max + 1 + i,
                                column_id=c['grid_column_name']))
    CustomColumn.objects.bulk_create(new)


def _reverse_append_to_all_columns(apps, _schema_editor):
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    CustomColumn.objects.filter(custom_columns_collection__name='All columns',
                                column_id__in=grid_names).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0131_dbnsfp_v4_pathogenicity_fields'),
        ('snpdb', '0175_settingsoverride_initially_show_zygosity_table'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
