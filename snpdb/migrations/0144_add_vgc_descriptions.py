# Generated by Django 4.2.15 on 2024-09-09 06:24

from django.db import migrations


def _add_vgc_descriptions(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")

    VGC_DESCRIPTIONS = {
        "alphamissense_rankscore": "AlphaMissense is an AI model developed by DeepMind that predicts pathogenicity of missense variants. See <a href='https://www.science.org/doi/10.1126/science.adg7492'>Paper</a>",
        "mavedb_urn": "DB identifier from <a href='https://www.mavedb.org/'>MaveDB</a> - a public repository for datasets from Multiplexed Assays of Variant Effect (MAVEs), such as those generated by deep mutational scanning or massively parallel reporter assay experiments.",
        "mavedb_score": "Score from <a href='https://www.mavedb.org/'>MaveDB</a> - a public repository for datasets from Multiplexed Assays of Variant Effect (MAVEs), such as those generated by deep mutational scanning or massively parallel reporter assay experiments.",
    }

    for pk, description in VGC_DESCRIPTIONS.items():
        VariantGridColumn.objects.filter(pk=pk).update(description=description)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0143_one_off_rename_clinvar_variantgridcolumns'),
    ]

    operations = [
        migrations.RunPython(_add_vgc_descriptions, reverse_code=migrations.RunPython.noop)
    ]
