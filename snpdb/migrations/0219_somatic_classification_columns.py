from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data

_NEW_COLUMNS = [
    {'grid_column_name': 'internally_classified_somatic',
     'variant_column': 'internally_classified_somatic',
     'annotation_level': 'D',
     'width': None,
     'label': 'Internal classification (Somatic)',
     'description': "Internal somatic clinical significances (AMP tiers) joined by '|'",
     'model_field': False,
     'queryset_field': False},
    {'grid_column_name': 'max_internal_somatic_classification',
     'variant_column': 'max_internal_somatic_classification',
     'annotation_level': 'D',
     'width': None,
     'label': 'Internal classification (Somatic Max)',
     'description': 'Highest internal somatic clinical significance (AMP tier)',
     'model_field': False,
     'queryset_field': False},
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
        ('snpdb', '0218_alter_vcflengthstats_variant_class'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
