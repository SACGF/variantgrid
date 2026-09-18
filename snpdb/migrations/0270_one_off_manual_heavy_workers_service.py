from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _somalier_enabled(_apps):
    return settings.SOMALIER.get("enabled")


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0269_all_variants_filter_gene_level_types"),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "Install and start the celeryd_heavy_workers service (somalier now runs there) - see https://github.com/SACGF/variantgrid/issues/1885",
        ], test=_somalier_enabled),
    ]
