from django.db import migrations
from django.db.models import Max

from library.django_utils.seed_data import insert_seed_data

_TRANSCRIPT = 'T'
_UTR_ANNOTATOR = 'https://github.com/Ensembl/UTRannotator'

_NEW_COLUMNS = [
    {'grid_column_name': 'utr_annotator_consequence',
     'variant_column': 'variantannotation__utr_annotator_consequence',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR uORF consequence",
     'description': ("UTRAnnotator effect of a 5' UTR variant on upstream ORFs: uAUG gained/lost, uORF stop "
                     f"codon gained/lost, or uORF frameshift. {_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_annotation',
     'variant_column': 'variantannotation__utr_annotator_annotation',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR uORF details",
     'description': ("UTRAnnotator detail of each uORF effect (&-joined), eg Kozak context and strength, "
                     "uORF type, distance to the CDS and Evidence (whether the uORF data file has it). "
                     f"{_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_existing_uorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_uorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing uORFs",
     'description': f"UTRAnnotator count of uORFs with a stop codon within the 5' UTR. {_UTR_ANNOTATOR}"},
    {'grid_column_name': 'utr_annotator_existing_out_of_frame_oorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_out_of_frame_oorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing out-of-frame oORFs",
     'description': ("UTRAnnotator count of out-of-frame ORFs in the 5' UTR that overlap the CDS. "
                     f"{_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_existing_inframe_oorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_inframe_oorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing in-frame oORFs",
     'description': f"UTRAnnotator count of in-frame ORFs in the 5' UTR that overlap the CDS. {_UTR_ANNOTATOR}"},
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
        ('annotation', '0004_utr_annotator'),
        ('snpdb', '0002_seed_data'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
