from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _check_has_unbackfilled_spliceai(apps):
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    any_ds_set = (
        Q(spliceai_pred_ds_ag__isnull=False)
        | Q(spliceai_pred_ds_al__isnull=False)
        | Q(spliceai_pred_ds_dg__isnull=False)
        | Q(spliceai_pred_ds_dl__isnull=False)
    )
    return VariantAnnotation.objects.filter(any_ds_set, spliceai_max_ds__isnull=True).exists()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0135_spliceai_fields"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_historical_spliceai_max_ds"]),
                        test=_check_has_unbackfilled_spliceai),
    ]
