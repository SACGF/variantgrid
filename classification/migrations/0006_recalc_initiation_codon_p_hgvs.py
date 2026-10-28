from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _test_initiation_codon_p_hgvs(apps):
    """ Start codon changes were stored as substitutions (ClinGen's p.Met1Val) or failed to compute (#1383) """
    ResolvedVariantInfo = apps.get_model("classification", "ResolvedVariantInfo")
    met1_change_q = Q(p_hgvs__regex=r"p\.\(Met1[A-Z_]")
    return ResolvedVariantInfo.objects.filter(met1_change_q | Q(p_hgvs_error__contains="AARefAlt")).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0005_remove_old_flags'),
        ('manual', '0001_initial'),
    ]

    operations = [
        ManualOperation.operation_manage(["recalc_resolved_variant_info_p_hgvs"], test=_test_initiation_codon_p_hgvs),
    ]
