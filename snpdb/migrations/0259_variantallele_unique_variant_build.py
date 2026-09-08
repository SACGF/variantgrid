"""
A variant can only ever have 1 Allele per build - see https://github.com/SACGF/variantgrid/issues/1361

The old ("variant", "genome_build", "allele") unique_together let a fresh Allele slip past
bulk_create(ignore_conflicts=True) every time populate_clingen_alleles_for_variants ran over a variant that
can never get a ClinGen Allele. Dedupe what leaked, then narrow the constraint so it can't happen again.
"""
import logging

from django.db import migrations
from django.db.models import Count, F

from manual.operations.manual_operations import ManualOperation


def _pick_keeper(variant_alleles, clingen_allele_ids, referenced_allele_ids):
    """ A ClinGen Allele is the real identity, then whatever classifications/tags already point at,
        then lowest pk - which is what Variant.allele has been returning all along """
    for va in variant_alleles:
        if va.allele_id in clingen_allele_ids:
            return va
    for va in variant_alleles:
        if va.allele_id in referenced_allele_ids:
            return va
    return variant_alleles[0]


def _dedupe_variant_alleles(apps, schema_editor):
    VariantAllele = apps.get_model("snpdb", "VariantAllele")
    AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")
    VariantTag = apps.get_model("analysis", "VariantTag")
    Classification = apps.get_model("classification", "Classification")
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    ClinicalContext = apps.get_model("classification", "ClinicalContext")
    ClinVarRecordCollection = apps.get_model("annotation", "ClinVarRecordCollection")

    dupe_keys = (VariantAllele.objects.values("variant_id", "genome_build_id")
                 .annotate(num_alleles=Count("pk")).filter(num_alleles__gt=1))

    unresolvable = []
    for key in dupe_keys:
        variant_alleles = list(VariantAllele.objects.filter(variant_id=key["variant_id"],
                                                            genome_build_id=key["genome_build_id"])
                               .select_related("allele").order_by("pk"))
        allele_ids = [va.allele_id for va in variant_alleles]
        clingen_allele_ids = {va.allele_id for va in variant_alleles if va.allele.clingen_allele_id}
        referenced_allele_ids = set(Classification.objects.filter(allele__in=allele_ids)
                                    .values_list("allele_id", flat=True))
        referenced_allele_ids.update(VariantTag.objects.filter(allele__in=allele_ids)
                                     .values_list("allele_id", flat=True))

        keeper_va = _pick_keeper(variant_alleles, clingen_allele_ids, referenced_allele_ids)
        keeper = keeper_va.allele
        for va in variant_alleles:
            if va.pk == keeper_va.pk:
                continue
            other = va.allele
            if other.clingen_allele_id and keeper.clingen_allele_id:
                # One variant can't register 2 ClinGen Allele IDs - stop rather than guess which is right
                unresolvable.append((key["variant_id"], key["genome_build_id"], keeper.pk, other.pk))
                continue

            if other.flag_collection_id and not keeper.flag_collection_id:
                keeper.flag_collection_id = other.flag_collection_id
                keeper.save(update_fields=["flag_collection_id"])

            num_tags = VariantTag.objects.filter(allele=other).update(allele=keeper)
            num_classifications = Classification.objects.filter(allele=other).update(allele=keeper)
            ImportedAlleleInfo.objects.filter(allele=other).update(allele=keeper)
            ClinVarRecordCollection.objects.filter(allele=other).update(allele=keeper)
            # AlleleLiftover is unique on (liftover, allele) - both alleles were often in the same run, and
            # the keeper's record of it is the one to keep
            keeper_liftover_ids = set(AlleleLiftover.objects.filter(allele=keeper)
                                      .values_list("liftover_id", flat=True))
            AlleleLiftover.objects.filter(allele=other).exclude(liftover_id__in=keeper_liftover_ids) \
                .update(allele=keeper)

            # ClinicalContext is unique on (allele, allele_origin_bucket, name) - move the ones that don't
            # collide and leave the rest where they are, exactly as Allele.merge() does. Never delete one -
            # discordance history hangs off it. 'classification_fix_allele_links' re-homes what's left.
            keeper_cc_keys = set(ClinicalContext.objects.filter(allele=keeper)
                                 .values_list("allele_origin_bucket", "name"))
            for cc in ClinicalContext.objects.filter(allele=other):
                if (cc.allele_origin_bucket, cc.name) not in keeper_cc_keys:
                    cc.allele = keeper
                    cc.save(update_fields=["allele"])

            # The surplus Allele's links in other builds come across too (mirrors Allele.merge)
            keeper_build_ids = set(VariantAllele.objects.filter(allele=keeper)
                                   .values_list("genome_build_id", flat=True))
            for other_va in VariantAllele.objects.filter(allele=other).exclude(pk=va.pk):
                if other_va.genome_build_id in keeper_build_ids:
                    AlleleLiftover.objects.filter(variant_allele=other_va).update(variant_allele=None)
                    other_va.delete()
                else:
                    other_va.allele = keeper
                    other_va.save(update_fields=["allele"])
                    keeper_build_ids.add(other_va.genome_build_id)

            AlleleLiftover.objects.filter(variant_allele=va).update(variant_allele=None)
            va.delete()
            logging.info("Variant %s (build %s): kept allele %s, merged %s (%d tags, %d classifications)",
                         key["variant_id"], key["genome_build_id"], keeper.pk, other.pk,
                         num_tags, num_classifications)

    if unresolvable:
        raise RuntimeError("VariantAllele duplicates where both Alleles have a ClinGenAllele need a human "
                           f"(variant, build, keeper allele, surplus allele): {unresolvable}")


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
