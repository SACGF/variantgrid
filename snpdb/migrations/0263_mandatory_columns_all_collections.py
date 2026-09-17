"""
Tags, tags_global and Sample became mandatory #1199 - put them back into any collection a user had
removed them from since 0246 made them lead every collection. Appended at the end so nothing else moves.

Reverse leaves them in place - which were missing isn't recorded anywhere.
"""
from django.db import migrations
from django.db.models import Max

_MANDATORY_COLUMNS = ["variant", "tags", "tags_global", "Sample"]


def _add_missing_mandatory_columns(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.all():
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        existing = set(cc_qs.values_list("column_id", flat=True))
        missing = [c for c in _MANDATORY_COLUMNS if c not in existing]
        if not missing:
            continue

        next_sort_order = (cc_qs.aggregate(max=Max("sort_order"))["max"] or 0) + 1
        for i, column_id in enumerate(missing):
            CustomColumn.objects.create(custom_columns_collection=ccc, column_id=column_id,
                                        sort_order=next_sort_order + i)
        # Historical models skip CustomColumn.save(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0262_duo_sibling"),
    ]

    operations = [
        migrations.RunPython(_add_missing_mandatory_columns, reverse_code=migrations.RunPython.noop),
    ]
