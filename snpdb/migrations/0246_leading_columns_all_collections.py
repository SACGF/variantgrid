"""
The special columns lead every collection, not just the system ones #686

Reverse leaves the order alone - a collection's previous arrangement isn't recorded anywhere.
"""
from django.db import migrations

# What a reader scans first, in order. Sample marks where the per-sample zygosity columns go
_LEADING_COLUMNS = ["variant", "classifications", "tags", "tags_global", "Sample"]


def _lead_with_special_columns(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.all():
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        ordered = list(cc_qs.order_by("sort_order").values_list("column_id", flat=True))
        for column_id in _LEADING_COLUMNS:
            if column_id not in ordered:
                CustomColumn.objects.create(custom_columns_collection=ccc, column_id=column_id,
                                            sort_order=len(ordered))
                ordered.append(column_id)

        new_order = _LEADING_COLUMNS + [c for c in ordered if c not in _LEADING_COLUMNS]
        for sort_order, column_id in enumerate(new_order):
            cc_qs.filter(column_id=column_id).update(sort_order=sort_order)

        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0245_classifications_column_description'),
    ]

    operations = [
        migrations.RunPython(_lead_with_special_columns, reverse_code=migrations.RunPython.noop),
    ]
