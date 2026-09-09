"""
annotsv_acmg_class is a Pathogenicity choice field now, so the grid, the CSV and the column summary
all show the class as a word - the description no longer needs the 1..5 code table.
"""
from django.db import migrations

_ANNOTSV_LINK = "<a href='https://lbgi.fr/AnnotSV/' target='_blank'>AnnotSV</a>"

_ACMG_CLASS = 'annotsv_acmg_class'

_DESCRIPTION = (
    f"{_ANNOTSV_LINK} ACMG-style ranking class for SVs: Benign, Likely benign, Uncertain, "
    f"Likely pathogenic or Pathogenic."
)

_DESCRIPTION_0180 = (
    f"{_ANNOTSV_LINK} ACMG-style ranking class for SVs (1=benign, 2=likely benign, 3=VUS, "
    f"4=likely pathogenic, 5=pathogenic)."
)


def _update_description(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk=_ACMG_CLASS).update(description=_DESCRIPTION)


def _restore_description(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk=_ACMG_CLASS).update(description=_DESCRIPTION_0180)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0259_variantallele_unique_variant_build'),
    ]

    operations = [
        migrations.RunPython(_update_description, reverse_code=_restore_description),
    ]
