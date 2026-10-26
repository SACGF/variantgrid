from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_cohort_genotypes(apps):
    return apps.get_model("snpdb", "CohortGenotypeCollection").objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0001_initial'),
        ('snpdb', '0005_test_variantgrid_site'),
    ]

    operations = [
        ManualOperation.operation_manage(["one_off_dedupe_cohort_genotypes"], test=_test_has_cohort_genotypes),
    ]
