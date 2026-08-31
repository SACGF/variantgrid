from django.db import migrations

# Now drawn inside the gnomAD AF cell. Stays in the catalogue and 'All columns', so it is still
# available from the column chooser and still exports
_MERGED_INTO_COMPOSITE = ["gnomad_filtered"]
_KEEP_EVERYTHING = "All columns"


def _update_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="gnomad_af").update(
        width=130,
        description="Allele Frequency (0-1) among all gnomAD genotypes (exome+genome), with the gnomAD "
                    "Pass/Fail link alongside it. Sorts by frequency")


def _update_system_collections(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
        if ccc.name != _KEEP_EVERYTHING:
            CustomColumn.objects.filter(custom_columns_collection=ccc,
                                        column_id__in=_MERGED_INTO_COMPOSITE).delete()
        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0225_classifications_column_and_representative_variant'),
    ]

    operations = [
        migrations.RunPython(_update_columns, reverse_code=migrations.RunPython.noop),
        migrations.RunPython(_update_system_collections, reverse_code=migrations.RunPython.noop),
    ]
