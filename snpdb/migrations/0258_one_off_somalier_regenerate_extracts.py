from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_somalier_extracts(apps):
    SomalierVCFExtract = apps.get_model("snpdb", "SomalierVCFExtract")
    return SomalierVCFExtract.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0257_one_off_fail_samples_of_failed_vcfs"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["somalier_existing_vcfs", "--clear"]),
                        note="Somalier now exports real allele depths (or GT alone) so it can QC the calls "
                             "itself - existing extracts, ancestry runs and relate pairs were built from the "
                             "old broken depths and need regenerating (#183)",
                        test=_has_somalier_extracts),
    ]
