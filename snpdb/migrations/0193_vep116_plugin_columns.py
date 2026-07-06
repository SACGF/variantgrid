from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data


_VARIANT = 'V'

_PROTVAR = 'https://www.ebi.ac.uk/ProtVar/'
_OPEN_TARGETS = 'https://platform.opentargets.org/'
_EVE = 'https://evemodel.org/'

_NEW_COLUMNS = [
    # ---- ProtVar ----
    {'grid_column_name': 'protvar_stability',
     'variant_column': 'variantannotation__protvar_stability',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'ProtVar ddG stability',
     'description': f'Impact on protein stability (ddG, kcal/mol; >2 likely destabilising). {_PROTVAR}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'protvar_pocket',
     'variant_column': 'variantannotation__protvar_pocket',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'ProtVar pocket',
     'description': ('Overlapping protein pocket (&-joined: id, score, mean pLDDT, energy, '
                     f'buriedness, radius of gyration, residues). {_PROTVAR}'),
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'protvar_int',
     'variant_column': 'variantannotation__protvar_int',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'ProtVar interface',
     'description': ('Overlapping protein-protein interface (&-joined: partner protein UniProt id, '
                     f'pDockQ docking score). {_PROTVAR}'),
     'model_field': True, 'queryset_field': True},

    # ---- Open Targets ----
    {'grid_column_name': 'open_targets_gwas_l2g_score',
     'variant_column': 'variantannotation__open_targets_gwas_l2g_score',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets GWAS L2G score',
     'description': f'Locus-to-gene (L2G) score predicting the causal gene at a GWAS locus. {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_gwas_gene_id',
     'variant_column': 'variantannotation__open_targets_gwas_gene_id',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets GWAS gene',
     'description': f'Ensembl gene id of the GWAS-implicated target. {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_gwas_diseases',
     'variant_column': 'variantannotation__open_targets_gwas_diseases',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets GWAS diseases',
     'description': f'Diseases associated at this GWAS locus. {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_study_type',
     'variant_column': 'variantannotation__open_targets_study_type',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets study type',
     'description': f'Open Targets study type (e.g. gwas, qtl). {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_study_id',
     'variant_column': 'variantannotation__open_targets_study_id',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets study id',
     'description': 'Open Targets study identifier - look up at https://platform.opentargets.org/study/<id>',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_variant_id',
     'variant_column': 'variantannotation__open_targets_variant_id',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets variant id',
     'description': 'Open Targets variant identifier - look up at https://platform.opentargets.org/variant/<id>',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_qtl_gene_id',
     'variant_column': 'variantannotation__open_targets_qtl_gene_id',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets QTL gene',
     'description': f'Ensembl gene id of the QTL-implicated target. {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'open_targets_qtl_biosample',
     'variant_column': 'variantannotation__open_targets_qtl_biosample',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets QTL biosample',
     'description': f'Biosample/tissue for the QTL association. {_OPEN_TARGETS}',
     'model_field': True, 'queryset_field': True},

    # ---- EVE / popEVE ----
    {'grid_column_name': 'eve_score',
     'variant_column': 'variantannotation__eve_score',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'EVE score',
     'description': f'EVE (evolutionary model of variant effect) score for missense variants. {_EVE}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'eve_class',
     'variant_column': 'variantannotation__eve_class',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'EVE class',
     'description': f'EVE classification (Benign / Uncertain / Pathogenic; 75% uncertain threshold). {_EVE}',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'popeve_score',
     'variant_column': 'variantannotation__popeve_score',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'popEVE score',
     'description': 'popEVE score (population-adjusted EVE) for missense variants.',
     'model_field': True, 'queryset_field': True},

    # ---- PromoterAI ----
    {'grid_column_name': 'promoter_ai_score',
     'variant_column': 'variantannotation__promoter_ai_score',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'PromoterAI score',
     'description': 'PromoterAI predicted impact of a promoter variant on gene expression (Illumina).',
     'model_field': True, 'queryset_field': True},
    {'grid_column_name': 'promoter_ai_tss_pos',
     'variant_column': 'variantannotation__promoter_ai_tss_pos',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'PromoterAI TSS position',
     'description': 'Transcription start site position used by PromoterAI.',
     'model_field': True, 'queryset_field': True},
]


def _new_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COLUMNS)])


def _reverse_new_columns(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
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
        ('annotation', '0156_vep116_plugin_fields'),
        ('snpdb', '0192_usersettingsoverride_loading_animations'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
