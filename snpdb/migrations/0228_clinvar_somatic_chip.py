from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data

_ALL_COLUMNS = "All columns"

# The oncogenicity text, the somatic twin of 'clinical_significance' - the ClinVar somatic chip's
# tooltip reads it, and it belongs in the catalogue beside the columns migration 0222 added
_ONCOGENIC_CLASSIFICATION_COLUMN = {
    'grid_column_name': 'clinvar_oncogenic_classification',
    'variant_column': 'clinvar__oncogenic_classification',
    'annotation_level': 'C',
    'width': None,
    'label': 'ClinVar Oncogenic Classification',
    'description': "ClinVar oncogenicity classification text, from the ONC field",
    'model_field': True,
    'queryset_field': True,
}


def _add_column(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", [_ONCOGENIC_CLASSIFICATION_COLUMN])])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="classifications").update(
        width=230,
        description='ClinVar and internal classifications, germline pair then somatic pair. Hover a '
                    'chip for the per-record summary, click it to scroll to the underlying column')

    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    all_columns = CustomColumnsCollection.objects.filter(name=_ALL_COLUMNS).first()
    if all_columns:
        sort_order_max = all_columns.customcolumn_set.aggregate(Max("sort_order"))["sort_order__max"] or 0
        CustomColumn.objects.get_or_create(custom_columns_collection=all_columns,
                                           column_id="clinvar_oncogenic_classification",
                                           defaults={"sort_order": sort_order_max + 1})

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


def _remove_column(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="clinvar_oncogenic_classification").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0227_composite_columns_spliceai_maxentscan_mastermind'),
    ]

    operations = [
        migrations.RunPython(_add_column, reverse_code=_remove_column),
    ]
