from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_variant_tags_without_allele(apps):
    return apps.get_model("analysis", "VariantTag").objects.filter(allele__isnull=True).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0003_seed_data'),
        ('manual', '0001_initial'),
    ]

    operations = [
        ManualOperation.operation_manage(["variant_tags", "link-alleles"],
                                         test=_test_has_variant_tags_without_allele),
    ]
