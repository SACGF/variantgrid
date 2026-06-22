from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data


_VARIANT = 'V'

_DESC_PARALLEL = ('&-separated parallel array (one element per contributing denovo-db record). '
                  'See https://denovo-db.gs.washington.edu/.')

_NEW_COLUMNS = [
    {'grid_column_name': 'denovo_db_studies',
     'variant_column': 'variantannotation__denovo_db_studies',
     'annotation_level': _VARIANT,
     'width': None,
     'label': 'denovo-db studies',
     'description': f'Contributing denovo-db studies. {_DESC_PARALLEL}',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'denovo_db_pubmed_ids',
     'variant_column': 'variantannotation__denovo_db_pubmed_ids',
     'annotation_level': _VARIANT,
     'width': None,
     'label': 'denovo-db PubMed IDs',
     'description': f'PubMed IDs of contributing denovo-db records. {_DESC_PARALLEL}',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'denovo_db_primary_phenotypes',
     'variant_column': 'variantannotation__denovo_db_primary_phenotypes',
     'annotation_level': _VARIANT,
     'width': None,
     'label': 'denovo-db primary phenotypes',
     'description': f'Primary phenotype per contributing denovo-db record. {_DESC_PARALLEL}',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'denovo_db_case_count',
     'variant_column': 'variantannotation__denovo_db_case_count',
     'annotation_level': _VARIANT,
     'width': None,
     'label': 'denovo-db case count',
     'description': 'Number of denovo-db records with a non-control PrimaryPhenotype.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'denovo_db_control_count',
     'variant_column': 'variantannotation__denovo_db_control_count',
     'annotation_level': _VARIANT,
     'width': None,
     'label': 'denovo-db control count',
     'description': 'Number of denovo-db records with PrimaryPhenotype=control.',
     'model_field': True,
     'queryset_field': True},
]

_NEW_COLUMN_VCF_INFO = [
    {'info_id': 'denovo_db_studies',              'column_id': 'denovo_db_studies',              'number': None, 'type': 'S', 'description': 'denovo-db contributing studies (&-separated)'},
    {'info_id': 'denovo_db_pubmed_ids',           'column_id': 'denovo_db_pubmed_ids',           'number': None, 'type': 'S', 'description': 'denovo-db PubMed IDs (&-separated)'},
    {'info_id': 'denovo_db_primary_phenotypes',   'column_id': 'denovo_db_primary_phenotypes',   'number': None, 'type': 'S', 'description': 'denovo-db primary phenotypes (&-separated)'},
    {'info_id': 'denovo_db_case_count',           'column_id': 'denovo_db_case_count',           'number': 1,    'type': 'I', 'description': 'denovo-db case count (non-control)'},
    {'info_id': 'denovo_db_control_count',        'column_id': 'denovo_db_control_count',        'number': 1,    'type': 'I', 'description': 'denovo-db control count'},
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
        ('annotation', '0133_denovo_db_fields'),
        ('snpdb', '0176_dbnsfp_v4_pathogenicity_columns'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
