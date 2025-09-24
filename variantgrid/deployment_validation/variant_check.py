from django.conf import settings

from snpdb.models import Variant


def check_symbolic_variants() -> dict:
    qs = Variant.objects.filter()
    fix = None
    if settings.VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE:
        qs = qs.filter(svlen__lt=0)
        if count := qs.count():
            fix = f"{settings.VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE=}, but there are {count} Variants with svlen<0"
    else:
        qs_del_pos_svlen = qs.filter(alt__seq='<DEL>', svlen__gt=0)
        if count := qs_del_pos_svlen.count():
            fix = f"{settings.VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE=}, but there are {count} <DEL> variants with svlen>0"

    data = {
        "variant_grid_columns": {
            "valid": not bool(fix),
            "fix": fix,
        }
    }
    return data
