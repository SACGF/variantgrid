from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_trios(apps):
    Trio = apps.get_model("snpdb", "Trio")
    return Trio.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0153_samplenode_copy_ratio_thresholds"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_recalculate_trio_karyomapping"]),
                        note="Recount trio karyomapping now reference calls are binned (#1754)",
                        test=_has_trios),
    ]
