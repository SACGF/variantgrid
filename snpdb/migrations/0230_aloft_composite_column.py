from django.db import migrations

_KEEP_EVERYTHING = "All columns"

_ALOFT_PRED = "aloft_pred"
_ALOFT_PARTNERS = ["aloft_high_confidence", "aloft_prob_tolerant", "aloft_prob_recessive",
                   "aloft_prob_dominant", "aloft_ensembl_transcript"]

_ALOFT_PRED_UPDATE = (
    125,
    "ALoFT",
    "<a href='https://www.nature.com/articles/s41467-017-00443-5'>ALoFT</a> classification of a "
    "loss of function variant - tolerant, recessive or dominant - with the probability behind it. "
    "A call ALoFT is not confident in (p &gt;= 0.05) is greyed. Hover for all three probabilities "
    "and the transcript it chose",
)


def _update_aloft_column(apps, _schema_editor):
    width, label, description = _ALOFT_PRED_UPDATE
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk=_ALOFT_PRED).update(width=width, label=label,
                                                            description=description)


def _update_system_collections(apps, _schema_editor):
    """ The partners leave the system collections and stay in the catalogue and 'All columns' """
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True).exclude(name=_KEEP_EVERYTHING):
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        ordered = list(cc_qs.order_by("sort_order").values_list("column_id", flat=True))
        group = [c for c in ordered if c == _ALOFT_PRED or c in _ALOFT_PARTNERS]
        if not group:
            continue
        position = ordered.index(group[0])
        ordered = [c for c in ordered if c not in _ALOFT_PARTNERS and c != _ALOFT_PRED]
        ordered.insert(position, _ALOFT_PRED)

        for column_id in set(ordered).difference(cc_qs.values_list("column_id", flat=True)):
            CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column_id=column_id,
                                               defaults={"sort_order": 0})
        cc_qs.exclude(column_id__in=ordered).delete()
        for sort_order, column_id in enumerate(ordered):
            cc_qs.filter(column_id=column_id).update(sort_order=sort_order)

        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0229_settingsoverride_variant_grid_two_line_rows'),
    ]

    operations = [
        migrations.RunPython(_update_aloft_column, reverse_code=migrations.RunPython.noop),
        migrations.RunPython(_update_system_collections, reverse_code=migrations.RunPython.noop),
    ]
