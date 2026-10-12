from django.conf import settings
from django.db import migrations
from django.db.models import Q
from django.db.models.functions import Length

from manual.operations.manual_operations import ManualOperation

_EXPLICIT_SYMBOLIC_ALTS = ["<DEL>", "<DUP>", "<INV>"]


def _test_non_canonical_symbolic_size(apps):
    """ VARIANT_SYMBOLIC_ALT_SIZE went from 1000 to 50 (#1358), and a symbolic variant could be stored just
        under the old threshold (#2109) """
    if not settings.VARIANT_SYMBOLIC_ALT_ENABLED:
        return False

    Sequence = apps.get_model("snpdb", "Sequence")
    Variant = apps.get_model("snpdb", "Variant")
    size = settings.VARIANT_SYMBOLIC_ALT_SIZE

    short_symbolic = Variant.objects.filter(alt__seq__in=_EXPLICIT_SYMBOLIC_ALTS, svlen__gt=-size, svlen__lt=size)
    if short_symbolic.exists():
        return True

    long_sequences = Sequence.objects.annotate(seq_length=Length("seq")).filter(seq_length__gt=size)
    long_explicit = Variant.objects.filter(Q(locus__ref__in=long_sequences) | Q(alt__in=long_sequences),
                                           svlen__isnull=True)
    return long_explicit.exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0278_remove_user_award_titles'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_fix_symbolic_variants"]),
                        note="Run with --dry-run first to see the counts",
                        test=_test_non_canonical_symbolic_size),
    ]
