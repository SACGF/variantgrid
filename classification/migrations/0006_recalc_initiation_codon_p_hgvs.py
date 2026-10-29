from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _test_p_hgvs_needs_recalc(apps):
    """ Start codon changes were stored as substitutions (ClinGen's p.Met1Val) or failed to compute (#1383).
        Also rows never backfilled - every calculation sets p_hgvs or p_hgvs_error """
    ResolvedVariantInfo = apps.get_model("classification", "ResolvedVariantInfo")
    met1_change_q = Q(p_hgvs__regex=r"p\.\(Met1[A-Z_]") | Q(p_hgvs_error__contains="AARefAlt")
    never_calculated_q = Q(transcript_version__isnull=False, p_hgvs__isnull=True, p_hgvs_error__isnull=True)
    return ResolvedVariantInfo.objects.filter(met1_change_q | never_calculated_q).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0005_remove_old_flags'),
        ('manual', '0001_initial'),
    ]

    operations = [
        # Computing p.HGVS needs cdot's protein accessions and the cached transcript sequences
        ManualOperation.operation_manage(["recalc_resolved_variant_info_p_hgvs"], test=_test_p_hgvs_needs_recalc,
                                         requires=["cdot-current", "transcript-sequences-loaded"]),
    ]
