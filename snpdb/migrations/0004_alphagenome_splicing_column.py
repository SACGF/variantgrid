from django.db import migrations
from django.db.models import Max

from library.django_utils.seed_data import insert_seed_data

_NEW_COLUMNS = [
    {'grid_column_name': 'alphagenome_splicing',
     'variant_column': 'variantannotation__alphagenome_splicing',
     'annotation_level': 'V',
     'label': 'AlphaGenome splicing',
     'description': ("AlphaGenome (Google DeepMind) merged splicing score, combining its splice site, splice site "
                     "usage and splice junction predictions. Higher is a bigger predicted effect on splicing. SNVs "
                     "only. Licensed for non-commercial use only, and not for diagnostic procedures. "
                     "https://deepmind.google.com/science/alphagenome")},
]


def _new_columns(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [("snpdb.VariantGridColumn", _NEW_COLUMNS)])


def _reverse_new_columns(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    apps.get_model("snpdb", "VariantGridColumn").objects.filter(grid_column_name__in=grid_names).delete()


def _append_to_all_columns(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    all_columns = CustomColumnsCollection.objects.filter(name='All columns', user__isnull=True).first()
    if all_columns is None:
        return
    sort_order_max = all_columns.customcolumn_set.aggregate(Max("sort_order"))["sort_order__max"] or 0
    CustomColumn.objects.bulk_create([
        CustomColumn(custom_columns_collection=all_columns, sort_order=sort_order_max + 1 + i,
                     column_id=c['grid_column_name'])
        for i, c in enumerate(_NEW_COLUMNS)
    ])
    # CustomColumn.save() does this for the grid caches, but bulk_create skips it
    all_columns.version_id += 1
    all_columns.save(update_fields=["version_id"])


def _reverse_append_to_all_columns(apps, _schema_editor):
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    CustomColumn.objects.filter(custom_columns_collection__name='All columns',
                                column_id__in=grid_names).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0007_alphagenome_splicing'),
        ('snpdb', '0003_utr_annotator_avi_splicevault_columns'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
