"""
Collapses VariantAllele duplicates so a variant has at most one Allele per build - see
https://github.com/SACGF/variantgrid/issues/1361

Owns the dedupe that snpdb migration 0259 runs before narrowing the unique_together, and stays runnable
afterwards for a database that grows a new duplicate. Entry points: `dedupe_variant_alleles(apps)` (the
migration hands over its historical app registry) and the `one_off_dedupe_variant_alleles` command.
"""
import logging

from django.apps import apps as django_apps
from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Count


class _DryRunRollback(Exception):
    pass


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


def _delete_variant_allele(apps, variant_allele):
    AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")
    AlleleLiftover.objects.filter(variant_allele=variant_allele).update(variant_allele=None)
    variant_allele.delete()


def _drop_liftover_artefacts(apps, variant_alleles, genome_build_id):
    """ Two registered ClinGen Alleles can't be one variant, so at least one link came from a bad liftover.
        A liftover that wrote 2 rows for the one run leaves the Allele with another variant in this build -
        that marks this link as the artefact, so drop it and leave the Allele for liftover to re-link. """
    VariantAllele = apps.get_model("snpdb", "VariantAllele")

    remaining = []
    for va in variant_alleles:
        also_linked_elsewhere = VariantAllele.objects.filter(allele_id=va.allele_id,
                                                             genome_build_id=genome_build_id) \
            .exclude(pk=va.pk).exists()
        if va.allele.clingen_allele_id and also_linked_elsewhere:
            logging.info("Variant %s (build %s): dropped liftover artefact link to allele %s",
                         va.variant_id, genome_build_id, va.allele_id)
            _delete_variant_allele(apps, va)
        else:
            remaining.append(va)
    return remaining


def dedupe_variant_alleles(apps) -> list:
    """ Returns the (variant, build, allele) links dropped without merging, for a human to look at """
    VariantAllele = apps.get_model("snpdb", "VariantAllele")
    AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")
    VariantTag = apps.get_model("analysis", "VariantTag")
    Classification = apps.get_model("classification", "Classification")
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    ClinicalContext = apps.get_model("classification", "ClinicalContext")
    ClinVarRecordCollection = apps.get_model("annotation", "ClinVarRecordCollection")

    dupe_keys = (VariantAllele.objects.values("variant_id", "genome_build_id")
                 .annotate(num_alleles=Count("pk")).filter(num_alleles__gt=1))

    dropped = []
    for key in dupe_keys:
        variant_alleles = list(VariantAllele.objects.filter(variant_id=key["variant_id"],
                                                            genome_build_id=key["genome_build_id"])
                               .select_related("allele").order_by("pk"))
        if len([va for va in variant_alleles if va.allele.clingen_allele_id]) > 1:
            variant_alleles = _drop_liftover_artefacts(apps, variant_alleles, key["genome_build_id"])
            if len(variant_alleles) < 2:
                continue

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
                # Both Alleles are real and neither link looks like an artefact. Drop the surplus link so the
                # unique constraint can go on - liftover can restore a link, a failed migration stops a deploy
                logging.warning("Variant %s (build %s): kept allele %s, dropped link to allele %s - both have "
                                "a ClinGenAllele, needs a human",
                                key["variant_id"], key["genome_build_id"], keeper.pk, other.pk)
                dropped.append((key["variant_id"], key["genome_build_id"], keeper.pk, other.pk))
                _delete_variant_allele(apps, va)
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
                    _delete_variant_allele(apps, other_va)
                else:
                    other_va.allele = keeper
                    other_va.save(update_fields=["allele"])
                    keeper_build_ids.add(other_va.genome_build_id)

            _delete_variant_allele(apps, va)
            logging.info("Variant %s (build %s): kept allele %s, merged %s (%d tags, %d classifications)",
                         key["variant_id"], key["genome_build_id"], keeper.pk, other.pk,
                         num_tags, num_classifications)
    return dropped


class Command(BaseCommand):
    """
        One Allele per variant per build (#1361) - collapse the duplicates that leaked in before the
        unique constraint, and any a bad liftover has added since
    """
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        try:
            with transaction.atomic():
                dropped = dedupe_variant_alleles(django_apps)
                for variant_id, genome_build_id, keeper_id, other_id in dropped:
                    print(f"Variant {variant_id} ({genome_build_id}): kept allele {keeper_id}, "
                          f"dropped link to allele {other_id} - both have a ClinGenAllele, needs a human")
                if dry_run:
                    raise _DryRunRollback()
        except _DryRunRollback:
            print("Dry run - rolled back")
