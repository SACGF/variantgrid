from django.db import migrations

from library.django_utils import bulk_insert_class_data

# Columns now drawn inside another cell. They stay in the catalogue and in 'All columns', so they are
# still available from the column chooser and still export
_REPLACED_BY_VARIANT_COLUMN = ["chrom", "position", "ref", "alt", "svlen", "hgvs_g"]
_MERGED_INTO_COMPOSITE = ["impact", "gnomad_popmax"]
_KEEP_EVERYTHING = "All columns"
# What a reader scans first, in order, straight after the mandatory Variant column
_LEADING_COLUMNS = ["variant", "classifications", "tags", "tags_global"]

_CLASSIFICATIONS_COLUMN = {
    'grid_column_name': 'classifications',
    'variant_column': 'classifications',
    'annotation_level': 'D',
    'width': 170,
    'label': 'Classifications',
    'description': 'ClinVar, internal germline and internal somatic classifications. Hover a chip for the '
                   'per-record summary, click it to scroll to the underlying column',
    'model_field': False,
    'queryset_field': False,
}


def _add_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", [_CLASSIFICATIONS_COLUMN])])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="variant").update(
        width=280,
        description="The variant named by c.HGVS, g.HGVS or coordinate - whichever is available. "
                    "Click the name for details, the row for transcript / HGVS / coordinates / IGV. "
                    "Sorts in genomic order")
    VariantGridColumn.objects.filter(pk="consequence").update(
        width=160,
        description="Sequence Ontology consequence, with a coloured dot for the impact "
                    "(MODIFIER/LOW/MODERATE/HIGH). Sorts by consequence")
    VariantGridColumn.objects.filter(pk="gnomad_popmax_af").update(
        width=110,
        description="Allele frequency in the outbred population with the highest AF (excludes Ashkenazi, "
                    "Finnish and other), with that population alongside it. Sorts by frequency")


def _update_system_collections(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        if ccc.name != _KEEP_EVERYTHING:
            cc_qs.filter(column_id__in=_REPLACED_BY_VARIANT_COLUMN + _MERGED_INTO_COMPOSITE).delete()
        CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column_id="classifications",
                                           defaults={"sort_order": 0})

        column_ids = list(cc_qs.order_by("sort_order").values_list("column_id", flat=True))
        leading = [c for c in _LEADING_COLUMNS if c in column_ids]
        for sort_order, column_id in enumerate(leading + [c for c in column_ids if c not in leading]):
            cc_qs.filter(column_id=column_id).update(sort_order=sort_order)

        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


def _remove_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="classifications").delete()  # cascades to CustomColumn


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0224_nodecountsettings_node_count_type'),
    ]

    operations = [
        migrations.RunPython(_add_columns, reverse_code=_remove_columns),
        migrations.RunPython(_update_system_collections, reverse_code=migrations.RunPython.noop),
    ]
