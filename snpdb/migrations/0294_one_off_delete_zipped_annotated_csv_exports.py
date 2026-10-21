from django.db import migrations

from manual.operations.manual_operations import ManualOperation

ANNOTATED_EXPORT_GENERATORS = ["export_cohort_to_downloadable_file", "export_sample_to_downloadable_file"]


def _has_zipped_annotated_csv_exports(apps):
    CachedGeneratedFile = apps.get_model("snpdb", "CachedGeneratedFile")
    return CachedGeneratedFile.objects.filter(generator__in=ANNOTATED_EXPORT_GENERATORS,
                                              filename__endswith=".csv.zip").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0293_one_off_somalier_hosted_build_rebuild_extracts"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_delete_zipped_annotated_csv_exports"]),
                        note="Delete zipped sample/cohort annotated CSV exports so pages offer the plain .csv (#2368)",
                        test=_has_zipped_annotated_csv_exports),
    ]
