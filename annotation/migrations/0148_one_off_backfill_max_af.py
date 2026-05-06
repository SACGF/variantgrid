import operator
from functools import reduce

from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _check_has_unbackfilled_max_af(apps):
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    for vav in VariantAnnotationVersion.objects.all():
        gnomad = vav.gnomad or ""
        try:
            major = int(gnomad.split(".", maxsplit=1)[0])
        except ValueError:
            continue
        fields = ["af_1kg", "af_uk10k", "gnomad_af", "gnomad_popmax_af"]
        if major >= 4:
            fields.append("gnomad_fafmax_faf95_max")
        all_set_q = reduce(operator.and_,
                           (Q(**{f"{f}__isnull": False}) for f in fields))
        if VariantAnnotation.objects.filter(version=vav, max_af__isnull=True).filter(all_set_q).exists():
            return True
    return False


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0147_variantannotation_max_af"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_historical_max_af"]),
                        test=_check_has_unbackfilled_max_af),
    ]
