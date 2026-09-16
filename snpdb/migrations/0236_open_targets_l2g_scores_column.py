from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data
from manual.operations.manual_operations import ManualOperation

_VARIANT = 'V'

_NEW_COLUMNS = [
    {'grid_column_name': 'open_targets_gwas_l2g_scores',
     'variant_column': 'variantannotation__open_targets_gwas_l2g_scores',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets GWAS L2G scores',
     'description': ('Per-record locus-to-gene scores, &-joined parallel to the other Open Targets '
                     'columns. https://platform.opentargets.org/'),
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


def _has_open_targets_annotation(apps):
    """ Only deployments that have annotated at columns_version 5 with an Open Targets data file
        have anything to backfill into the new per-record score column """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(columns_version__gte=5).exclude(open_targets__isnull=True) \
                                           .exclude(open_targets="").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0179_open_targets_gwas_l2g_scores'),
        ('snpdb', '0235_quad_proband_sex'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
        ManualOperation(task_id=ManualOperation.task_id_manage(["annotation_backfill_columns"]),
                        note="Backfill open_targets_gwas_l2g_scores for GRCh38 columns_version 5 annotation (#1822)",
                        test=_has_open_targets_annotation),
    ]
