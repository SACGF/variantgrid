from django.db import migrations
from django.db.models import Max

from library.django_utils import bulk_insert_class_data
from library.django_utils.composite_columns import collapse_into_composite

_VARIANT = 'V'
_COMPOSITE = 'open_targets'

_NEW_COLUMNS = [
    {'grid_column_name': 'open_targets_is_lead',
     'variant_column': 'variantannotation__open_targets_is_lead',
     'annotation_level': _VARIANT, 'width': None,
     'label': 'Open Targets lead variant',
     'description': ("Whether the variant leads each credible set it belongs to, &-joined parallel "
                     "to the other Open Targets columns. https://platform.opentargets.org/"),
     'model_field': True, 'queryset_field': True},
]


def _new_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COLUMNS)])


def _reverse_new_columns(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    apps.get_model("snpdb", "VariantGridColumn").objects.filter(grid_column_name__in=grid_names).delete()


def _add_composite_members(apps, _schema_editor):
    """ Detail only, like the rest of the Open Targets record fields - it reads on hover beside the
        study it qualifies rather than as a column of its own """
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    member_qs = CompositeColumnMember.objects.filter(composite_id=_COMPOSITE)
    sort_order = (member_qs.aggregate(Max("sort_order"))["sort_order__max"] or 0) + 1
    for i, c in enumerate(_NEW_COLUMNS):
        CompositeColumnMember.objects.create(composite_id=_COMPOSITE, column_id=c['grid_column_name'],
                                             sort_order=sort_order + i, in_sort_menu=False)
    collapse_into_composite(apps, _COMPOSITE)


def _reverse_composite_members(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    apps.get_model("snpdb", "CompositeColumnMember").objects.filter(composite_id=_COMPOSITE,
                                                                    column_id__in=grid_names).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0180_open_targets_is_lead'),
        ('snpdb', '0240_open_targets_l2g_scores_backfill_note'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_add_composite_members, reverse_code=_reverse_composite_members),
    ]
