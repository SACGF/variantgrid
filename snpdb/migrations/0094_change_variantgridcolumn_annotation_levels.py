# Generated by Django 4.1.3 on 2023-02-20 04:39

from django.db import migrations


def _change_variantgridcolumn_annotation_levels(apps, _schema_editor):
    variant_columns = [
        "variantannotation__dbsnp_rs_id",
        "variantannotation__phastcons_100_way_vertebrate",
        "variantannotation__phastcons_30_way_mammalian",
        "variantannotation__phastcons_46_way_mammalian",
        "variantannotation__phylop_100_way_vertebrate",
        "variantannotation__phylop_30_way_mammalian",
        "variantannotation__phylop_46_way_mammalian",
        "variantannotation__pubmed",
    ]
    variant_level = 'V'
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    vgc_qs = VariantGridColumn.objects.filter(variant_column__in=variant_columns)
    vgc_qs.update(annotation_level=variant_level)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0093_change_omim_id_label'),
    ]

    operations = [
        migrations.RunPython(_change_variantgridcolumn_annotation_levels),
    ]
