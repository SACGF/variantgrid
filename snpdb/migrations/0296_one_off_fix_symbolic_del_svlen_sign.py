from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_wrong_sign_svlen(apps):
    """ A <DEL> imported with END but no SVLEN was stored with a positive svlen before #1340 - a duplicate
        of the negative form. Mirrors one_off_fix_symbolic_variants.wrong_sign_svlen_q """
    if not settings.VARIANT_SYMBOLIC_ALT_ENABLED:
        return False

    Variant = apps.get_model("snpdb", "Variant")
    if settings.VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE:
        qs = Variant.objects.filter(svlen__lt=0)
    else:
        qs = Variant.objects.filter(alt__seq="<DEL>", svlen__gt=0)
    return qs.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0295_alter_allelemergelog_options_and_more'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_fix_symbolic_variants"]),
                        note="Run with --dry-run first to see the counts",
                        test=_test_wrong_sign_svlen),
    ]
