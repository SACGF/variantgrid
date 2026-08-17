from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation

_AFFECTED_CATEGORIES = ["liftover", "builds", "general"]


def _test_has_g_hgvs_only_validation_tags(apps):
    """ g.HGVS only submissions were validated as though a c.HGVS was expected - only worth running where
        such records exist and carry tags from one of the categories that changed (#1063) """
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    return ImportedAlleleInfo.objects.filter(
        Q(imported_c_hgvs__isnull=True) | Q(imported_c_hgvs=""),
        Q(imported_transcript__isnull=True) | Q(imported_transcript=""),
        latest_validation__validation_tags__has_any_keys=_AFFECTED_CATEGORIES,
    ).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0168_early_stop_codon_ekey'),
    ]

    operations = [
        ManualOperation.operation_manage(
            ["revalidate_g_hgvs_only_allele_infos"],
            note="Re-validate g.HGVS only classifications so they're no longer excluded from exports (#1063)",
            test=_test_has_g_hgvs_only_validation_tags),
    ]
