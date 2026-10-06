from django.db import migrations
from django.db.models import F


def _tumor_fraction_percent_to_fraction(apps, _schema_editor):
    """ sapath#457 - DRAGEN 2.6 writes the CombinedVariantOutput tumour fraction as a percent (55), which
        imports before the fix stored as is. A fraction is never above 1 """
    DragenTSO500CombinedVariantOutput = apps.get_model("seqauto", "DragenTSO500CombinedVariantOutput")
    DragenTSO500CombinedVariantOutput.objects.filter(tumor_fraction__gt=1).update(tumor_fraction=F("tumor_fraction") / 100)


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0050_seqauto_api_write_group"),
    ]

    operations = [
        migrations.RunPython(_tumor_fraction_percent_to_fraction, migrations.RunPython.noop),
    ]
