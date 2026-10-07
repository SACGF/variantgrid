from django.conf import settings
from django.db import migrations
from django.db.models import Min, Q

from manual.operations.manual_operations import ManualOperation


def _has_svs_too_long_for_vep(apps):
    """ SVs over VEP's size cap were left without a VariantAnnotation row - #2104 """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    Variant = apps.get_model("snpdb", "Variant")

    if not VariantAnnotationVersion.objects.exists():
        return False
    sv_max_size = VariantAnnotationVersion.objects.aggregate(Min("sv_max_size"))["sv_max_size__min"]
    sv_max_size = sv_max_size or settings.ANNOTATION_VEP_SV_MAX_SIZE
    if not sv_max_size:
        return False
    return Variant.objects.filter(Q(svlen__gt=sv_max_size) | Q(svlen__lt=-sv_max_size)).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0187_delete_version_diffs'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_annotation_vep_too_long"]),
                        note="Write vep_skipped_reason=TOO_LONG rows for SVs too long for VEP (#2104)",
                        test=_has_svs_too_long_for_vep),
    ]
