"""
Variant tag maintenance.

Tag names are their primary key, so deployments accumulate case variants of the same tag (eg 'artefact' and
'Artefact') that split reporting and node filters. Repeats of one analysis tagging can no longer be written
(@see VariantTag.Meta - its migration deleted the ones made before it); global taggings can still repeat,
and delete-duplicates clears those.

A tag's allele is set by variant_tag_created_task, so a tag whose task never finished (worker killed, broker outage,
a ClinGen crash - #2354) is left with none, and doesn't show on the same variant in other builds. link-alleles makes
the missing VariantAlleles in bulk (one ClinGen call per batch), links the tags in SQL and lifts the alleles over.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging
from contextlib import contextmanager

from django.conf import settings
from django.core.management import BaseCommand
from django.db.models import Exists, Min, OuterRef, QuerySet, Subquery
from django.db.models.signals import post_delete

from analysis.models import VariantTag
from analysis.signals.signal_handlers import variant_tag_delete
from eventlog.models import create_event
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import (
    _create_variant_allele_with_new_allele,
    populate_clingen_alleles_for_variants,
)
from snpdb.clingen_allele_api import (
    ClinGenAlleleRegistryUnavailableException,
    ClinGenAlleleServerException,
)
from snpdb.enums import ImportSource
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Allele, GenomeBuild, Tag, Variant, VariantAllele
from snpdb.tag_operations import get_case_collision_groups, get_tag_usage, merge_tag

MERGE_CASE_COLLISIONS = "merge-case-collisions"
DELETE_DUPLICATES = "delete-duplicates"
LINK_ALLELES = "link-alleles"

# What makes one variant tag a repeat of another - the analysis tagging identity (@see VariantTag.Meta).
# A different user re-tagging is agreement, and a different sample is a different person's to-do, so both stay
DUPLICATE_FIELDS = ["variant_id", "tag_id", "analysis_id", "user_id", "sample_id"]


@contextmanager
def _variant_tag_delete_signal_disconnected():
    """ Every duplicate we delete leaves an identical row behind, so the analyses they belong to see no
        change - skip the per-row 'set tag nodes dirty' + celery task the delete signal would fire.
        Also what lets Django delete them as one statement rather than a row at a time """
    post_delete.disconnect(variant_tag_delete, sender=VariantTag)
    try:
        yield
    finally:
        post_delete.connect(variant_tag_delete, sender=VariantTag)


def get_duplicate_variant_tags_qs() -> QuerySet[VariantTag]:
    """ Repeats of a variant tag, keeping the earliest so re-tagging history keeps its original event """
    earliest_qs = VariantTag.objects.values(*DUPLICATE_FIELDS).annotate(earliest=Min("pk")).values("earliest")
    return VariantTag.objects.exclude(pk__in=earliest_qs)


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="subcommand", required=True)
        merge_parser = subparsers.add_parser(
            MERGE_CASE_COLLISIONS, help="Merge every set of tags that differ only by case, keeping the most used")
        duplicates_parser = subparsers.add_parser(
            DELETE_DUPLICATES,
            help="Delete variant tags repeating the same variant/tag/analysis/user/sample")
        link_alleles_parser = subparsers.add_parser(
            LINK_ALLELES, help="Link variant tags with no allele to their variant's allele, then lift them over")
        for subparser in [merge_parser, duplicates_parser, link_alleles_parser]:
            subparser.add_argument('--dry-run', action='store_true',
                                   help="Report what would change without doing it")

    def handle(self, *args, **options):
        sub_commands = {
            MERGE_CASE_COLLISIONS: self._merge_case_collisions,
            DELETE_DUPLICATES: self._delete_duplicates,
            LINK_ALLELES: self._link_alleles,
        }
        sub_commands[options["subcommand"]](options["dry_run"])

    @staticmethod
    def _merge_case_collisions(dry_run: bool):
        for tag_ids in get_case_collision_groups():
            tags = sorted((Tag.objects.get(pk=t) for t in tag_ids),
                          key=lambda t: get_tag_usage(t).total, reverse=True)
            surviving_tag = tags[0]
            for dying_tag in tags[1:]:
                if dry_run:
                    logging.info("Would merge '%s' (%s) into '%s'",
                                 dying_tag, get_tag_usage(dying_tag).description(), surviving_tag)
                else:
                    merge_tag(dying_tag, surviving_tag, admin_bot())

    @staticmethod
    def _delete_duplicates(dry_run: bool):
        qs = get_duplicate_variant_tags_qs()
        if dry_run:
            logging.info("Would delete %d duplicate variant tags", qs.count())
        else:
            with _variant_tag_delete_signal_disconnected():
                num_deleted, _ = qs.delete()
            # No single tag to hang this off, so it's an event rather than a tag operation
            create_event(admin_bot(), name="variant_tags_delete_duplicates",
                         details=f"Deleted {num_deleted} duplicate variant tags")
            logging.info("Deleted %d duplicate variant tags", num_deleted)

    @staticmethod
    def _link_alleles(dry_run: bool):
        tags_qs = VariantTag.objects.filter(allele__isnull=True)
        for genome_build in GenomeBuild.objects.filter(pk__in=tags_qs.values("genome_build")):
            build_tags_qs = tags_qs.filter(genome_build=genome_build)
            if dry_run:
                logging.info("%s: would link %d variant tags with no allele", genome_build, build_tags_qs.count())
                continue

            # Held, as linking empties build_tags_qs. Few tags miss their task, so this is small
            variant_ids = list(build_tags_qs.values_list("variant_id", flat=True).distinct())
            variants_qs = Variant.objects.filter(pk__in=variant_ids)
            _create_missing_variant_alleles(genome_build, variants_qs)

            variant_allele_qs = VariantAllele.objects.filter(genome_build=genome_build, variant=OuterRef("variant"))
            num_linked = build_tags_qs.filter(Exists(variant_allele_qs)).update(
                allele=Subquery(variant_allele_qs.values("allele")))
            logging.info("%s: linked %d variant tags to their allele, %d still have none",
                         genome_build, num_linked, build_tags_qs.count())

            build_variant_allele_qs = VariantAllele.objects.filter(genome_build=genome_build, variant__in=variants_qs)
            allele_qs = Allele.objects.filter(pk__in=build_variant_allele_qs.values("allele"))
            create_liftover_pipelines(admin_bot(), allele_qs, ImportSource.COMMAND_LINE, genome_build)


def _create_missing_variant_alleles(genome_build: GenomeBuild, variants_qs: QuerySet[Variant]):
    """ As the tag task does for one variant: when ClinGen fails (or isn't configured) the variant still gets an
        Allele, with a server error stored so ClinGen is retried the next time the variant is looked up """
    clingen_error = None
    if settings.CLINGEN_ALLELE_REGISTRY_LOGIN:
        try:
            populate_clingen_alleles_for_variants(genome_build, variants_qs)
        except (ClinGenAlleleServerException, ClinGenAlleleRegistryUnavailableException) as e:
            logging.error("%s: ClinGen Allele Registry failed, making Alleles without ClinGen: %s", genome_build, e)
            clingen_error = e.get_fake_api_response()

    variant_allele_qs = VariantAllele.objects.filter(genome_build=genome_build, variant=OuterRef("pk"))
    for variant in variants_qs.exclude(Exists(variant_allele_qs)):
        _create_variant_allele_with_new_allele(variant, genome_build, clingen_error=clingen_error)
