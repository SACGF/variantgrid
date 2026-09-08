"""
The last four standalone AnnotSV columns fold into cells - a new gene effect composite (exons
spanned, frameshift, regulatory element genes) and the pathogenic SV overlaps on the hover of the
existing AnnotSV ACMG cell.

Reverse is lossy in the same way as snpdb/0238: the composite row goes and a collection keeps only
whatever members it still had.
"""
from django.db import migrations

from library.django_utils import bulk_insert_class_data
from library.django_utils.composite_columns import collapse_into_composite

_VARIANT = 'V'

_GENE_EFFECT = 'annotsv_gene_effect'
_ACMG = 'annotsv_acmg'

_NEW_COMPOSITE_COLUMNS = [
    {'grid_column_name': _GENE_EFFECT, 'variant_column': _GENE_EFFECT, 'annotation_level': _VARIANT,
     'width': None, 'label': 'AnnotSV gene effect',
     'description': 'Number of exons the structural variant fully spans, with whether it introduces a '
                    'frameshift and the regulatory element genes it overlaps on hover. Sorts by exons '
                    'spanned',
     'model_field': False, 'queryset_field': False},
]

_ACMG_DESCRIPTION = (
    'AnnotSV ACMG class for the structural variant, with the score, the criteria it applied and the '
    'overlapping pathogenic SVs on hover'
)

_ACMG_DESCRIPTION_0238 = (
    'AnnotSV ACMG class for the structural variant, with the score and the criteria it applied on hover'
)

# (composite, column, sort_order, in_sort_menu)
_MEMBERS = [
    (_GENE_EFFECT, 'annotsv_exons_spanned', 0, True),
    (_GENE_EFFECT, 'annotsv_frameshift', 1, True),
    (_GENE_EFFECT, 'annotsv_re_gene', 2, False),
    (_ACMG, 'annotsv_pathogenic_overlaps', 3, False),
]


def _add_composite_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COMPOSITE_COLUMNS)])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk=_ACMG).update(description=_ACMG_DESCRIPTION)


def _add_members(apps, _schema_editor):
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    CompositeColumnMember.objects.bulk_create([
        CompositeColumnMember(composite_id=composite, column_id=column_id, sort_order=sort_order,
                              in_sort_menu=in_sort_menu)
        for composite, column_id, sort_order, in_sort_menu in _MEMBERS
    ])


def _collapse_collections(apps, _schema_editor):
    collapse_into_composite(apps, _GENE_EFFECT)
    collapse_into_composite(apps, _ACMG)


def _remove_composite_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    # Cascades to the members and to the CustomColumns showing the composite
    VariantGridColumn.objects.filter(pk=_GENE_EFFECT).delete()
    CompositeColumnMember.objects.filter(composite_id=_ACMG,
                                         column_id='annotsv_pathogenic_overlaps').delete()
    VariantGridColumn.objects.filter(pk=_ACMG).update(description=_ACMG_DESCRIPTION_0238)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0248_merge_20260904_1656'),
    ]

    operations = [
        migrations.RunPython(_add_composite_columns, reverse_code=_remove_composite_columns),
        migrations.RunPython(_add_members, reverse_code=migrations.RunPython.noop),
        migrations.RunPython(_collapse_collections, reverse_code=migrations.RunPython.noop),
    ]
