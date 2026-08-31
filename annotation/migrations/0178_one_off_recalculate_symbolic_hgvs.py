from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation

PLACEHOLDER_MESSAGES = [
    "HGVS not calculated due to length",  # AbstractVariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
    "Error creating HGVS",  # AbstractVariantAnnotation.SV_HGVS_ERROR_MESSAGE
]
RANGED_SYMBOLIC_ALTS = ["<DEL>", "<DUP>", "<INV>"]


def _has_symbolic_hgvs_placeholders(apps):
    """ Symbolic DEL/DUP/INV only became expressible from coordinates alone in #1571 - before that
        both messages went into the annotation where the HGVS string should be """
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    VariantTranscriptAnnotation = apps.get_model("annotation", "VariantTranscriptAnnotation")

    for klass, placeholder_q in [
        (VariantAnnotation, Q(hgvs_g__in=PLACEHOLDER_MESSAGES) | Q(hgvs_c__in=PLACEHOLDER_MESSAGES)),
        (VariantTranscriptAnnotation, Q(hgvs_c__in=PLACEHOLDER_MESSAGES)),
    ]:
        if klass.objects.filter(placeholder_q, variant__alt__seq__in=RANGED_SYMBOLIC_ALTS).exists():
            return True
    return False


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0177_one_off_backfill_clinvar_somatic"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_recalculate_symbolic_hgvs"]),
                        note="Recalculate g./c.HGVS for symbolic DEL/DUP/INV annotation left holding "
                             "'HGVS not calculated due to length' or 'Error creating HGVS' (#1571)",
                        test=_has_symbolic_hgvs_placeholders),
    ]
