# Generated by Django 4.2.11 on 2024-08-06 11:08

from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _get_iai_qs(apps):
    Variant = apps.get_model("snpdb", "Variant")
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    bad_norm_qs = Variant.objects.filter(svlen__isnull=False, modifiedimportedvariant__isnull=False)
    return ImportedAlleleInfo.objects.filter(Q(grch37__variant__in=bad_norm_qs) | Q(grch38__variant__in=bad_norm_qs))


def _needs_hard_rematch(apps):
    qs = _get_iai_qs(apps)
    return qs.exists()


def _get_note(apps):
    qs = _get_iai_qs(apps)
    instructions = None
    if imported_allele_ids := list(qs.values_list("pk", flat=True)):
        iais = ", ".join((str(pk) for pk in imported_allele_ids))
        instructions = f"Go to /admin/classification/importedalleleinfo/ and select {iais}"
    return instructions


class Migration(migrations.Migration):
    dependencies = [
        ("classification", "0152_allele_origin_confirmed_ekey"),
        ("snpdb", "0141_one_off_fix_variant_end2"),
    ]

    operations = [
        ManualOperation.operation_other(args="AFTER fix_variant_end - Rematch hard imported alleles to badly normalized symbolic variants. ",
                                        note=_get_note, test=_needs_hard_rematch),
    ]