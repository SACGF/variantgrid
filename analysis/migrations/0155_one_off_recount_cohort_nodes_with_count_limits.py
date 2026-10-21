from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _has_cohort_nodes_with_count_limits(apps):
    CohortNode = apps.get_model("analysis", "CohortNode")
    any_limit_q = Q()
    for count_name in ["het_or_hom", "ref", "het", "hom"]:  # AbstractZygosityCountNode.COUNT_COLUMNS
        any_limit_q |= Q(**{f"min_{count_name}_count__gt": 0}) | Q(**{f"max_{count_name}_count__isnull": False})
    return CohortNode.objects.filter(any_limit_q, accordion_panel=0, cohort__isnull=False).exists()  # 0 = COUNT


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0154_one_off_recalculate_trio_karyomapping"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_recount_cohort_nodes_with_count_limits"]),
                        note="Reload cohort nodes with zygosity count limits, whose counts came from the "
                             "cohort's unfiltered stats (#2370)",
                        test=_has_cohort_nodes_with_count_limits),
    ]
