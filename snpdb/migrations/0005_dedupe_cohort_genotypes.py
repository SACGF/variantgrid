from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_cohort_genotypes(apps):
    return apps.get_model("snpdb", "CohortGenotypeCollection").objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0001_initial'),
        ('snpdb', '0004_alphagenome_splicing_column'),
    ]

    operations = [
        ManualOperation.operation_manage(["dedupe_cohort_genotypes"], test=_test_has_cohort_genotypes),
    ]
