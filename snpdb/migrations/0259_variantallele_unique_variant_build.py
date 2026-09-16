"""
A variant can only ever have 1 Allele per build - see https://github.com/SACGF/variantgrid/issues/1361

The old ("variant", "genome_build", "allele") unique_together let a fresh Allele slip past
bulk_create(ignore_conflicts=True) every time populate_clingen_alleles_for_variants ran over a variant that
can never get a ClinGen Allele. Dedupe what leaked, then narrow the constraint so it can't happen again.

The dedupe itself lives in the 'one_off_dedupe_variant_alleles' command so it can be run and re-run against
a database by hand.
"""
from django.db import migrations
from django.db.models import F

from manual.operations.manual_operations import ManualOperation
from snpdb.management.commands.one_off_dedupe_variant_alleles import dedupe_variant_alleles


def _dedupe_variant_alleles(apps, schema_editor):
    dedupe_variant_alleles(apps)


def _has_classifications_on_wrong_allele(apps) -> bool:
    Classification = apps.get_model("classification", "Classification")
    return Classification.objects.filter(allele__isnull=False, clinical_context__isnull=False) \
        .exclude(clinical_context__allele=F("allele")).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0258_one_off_somalier_regenerate_extracts"),
        ("analysis", "0141_varianttag_varianttag_one_per_sample_in_analysis"),
        ("classification", "0177_evidence_key_copy_scope_and_allele_origin_values"),
        ("annotation", "0180_open_targets_is_lead"),
    ]

    operations = [
        migrations.RunPython(_dedupe_variant_alleles, migrations.RunPython.noop),
        # The dedupe's deletes leave deferred FK triggers pending, and Postgres won't ALTER a table that has
        # them - fire them now so the constraint swap below can run in the same transaction
        migrations.RunSQL("SET CONSTRAINTS ALL IMMEDIATE", migrations.RunSQL.noop),
        migrations.AlterUniqueTogether(
            name="variantallele",
            unique_together={("variant", "genome_build")},
        ),
        ManualOperation(task_id=ManualOperation.task_id_manage(["classification_fix_allele_links"]),
                        note="Re-home clinical contexts and groupings of classifications moved off "
                             "duplicate alleles (#1361)",
                        test=_has_classifications_on_wrong_allele),
    ]
