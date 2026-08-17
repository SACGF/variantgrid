from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data

_TRANSCRIPT = 'T'

# sacgf/variantgrid#579 - predicted stop codon position + PTC-aware NMD
_NEW_COLUMNS = [
    {'grid_column_name': 'ptc_distance_codons',
     'variant_column': 'variantannotation__ptc_distance_codons',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'Stop codon distance',
     'description': 'Codons from the first changed residue to the frameshift\'s new premature termination codon (PTC), read from the HGVS p. "fsTer&lt;n&gt;". 1 means the changed codon is itself the stop.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'ptc_last_junction_distance',
     'variant_column': 'variantannotation__ptc_last_junction_distance',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'PTC to last junction',
     'description': 'Nucleotides from the premature termination codon to the last exon-exon junction. Negative means the PTC lies in the final exon. Drives the 50nt NMD rule.',
     'model_field': True,
     'queryset_field': True},
    {'grid_column_name': 'nmd_escape_status',
     'variant_column': 'variantannotation__nmd_escape_status',
     'annotation_level': _TRANSCRIPT,
     'width': None,
     'label': 'NMD escape',
     'description': 'Nonsense mediated decay prediction anchored on the premature termination codon rather than the variant. Escapes NMD if any of: 1. PTC is in the last exon 2. PTC is within 50 bases of the last exon-exon junction 3. PTC falls in the first 100 coding bases 4. Transcript has only 1 exon (no introns).',
     'model_field': True,
     'queryset_field': True},
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
        ('annotation', '0162_ptc_annotation'),
        ('snpdb', '0199_remove_samplestatspassingfilter_code_version_and_more'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
