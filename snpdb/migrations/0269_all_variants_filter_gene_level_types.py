from django.db import migrations


def _split_fusion_type(apps, _schema_editor):
    """ The 'fusion' variant type used to select every gene-level variant - keep saved filters' results """
    AllVariantsFilter = apps.get_model("snpdb", "AllVariantsFilter")
    for all_variants_filter in AllVariantsFilter.objects.filter(filters__variant_types__contains=["fusion"]):
        variant_types = all_variants_filter.filters["variant_types"]
        for variant_type in ["copy_number", "splice"]:
            if variant_type not in variant_types:
                variant_types.append(variant_type)
        all_variants_filter.save()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0268_cohort_phenotype"),
    ]

    operations = [
        migrations.RunPython(_split_fusion_type, migrations.RunPython.noop),
    ]
