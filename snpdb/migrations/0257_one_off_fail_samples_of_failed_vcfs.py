from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_open_samples_on_failed_vcfs(apps):
    Sample = apps.get_model("snpdb", "Sample")
    return Sample.objects.filter(vcf__import_status="E").exclude(import_status__in=["E", "M", "D"]).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0256_cohort_genotype_stats_fusions_count"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_fail_samples_of_failed_vcfs"]),
                        note="Mark samples of failed VCFs as failed - they used to be released independently",
                        test=_has_open_samples_on_failed_vcfs),
    ]
