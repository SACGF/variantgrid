"""
Variant tag maintenance.

Tag names are their primary key, so deployments accumulate case variants of the same tag (eg 'artefact' and
'Artefact') that split reporting and node filters. Merging them can leave the same user tagging the same
variant with the same tag in the same analysis twice, which is also worth cleaning up on its own.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging
from contextlib import contextmanager

from django.core.management import BaseCommand
from django.db.models import Min, QuerySet
from django.db.models.signals import post_delete

from analysis.models import VariantTag
from analysis.signals.signal_handlers import variant_tag_delete
from eventlog.models import create_event
from library.guardian_utils import admin_bot
from snpdb.models import Tag
from snpdb.tag_operations import get_case_collision_groups, get_tag_usage, merge_tag, set_tag_nodes_dirty

MERGE_CASE_COLLISIONS = "merge-case-collisions"
DELETE_DUPLICATES = "delete-duplicates"

# What makes one variant tag a repeat of another. A different user re-tagging is agreement, so is kept
DUPLICATE_FIELDS = ["variant_id", "tag_id", "analysis_id", "user_id"]


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
    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="subcommand", required=True)
        merge_parser = subparsers.add_parser(
            MERGE_CASE_COLLISIONS, help="Merge every set of tags that differ only by case, keeping the most used")
        duplicates_parser = subparsers.add_parser(
            DELETE_DUPLICATES, help="Delete variant tags repeating the same variant/tag/analysis/user")
        for subparser in [merge_parser, duplicates_parser]:
            subparser.add_argument('--dry-run', action='store_true',
                                   help="Report what would change without doing it")

    def handle(self, *args, **options):
        sub_commands = {
            MERGE_CASE_COLLISIONS: self._merge_case_collisions,
            DELETE_DUPLICATES: self._delete_duplicates,
        }
        sub_commands[options["subcommand"]](options["dry_run"])

    @staticmethod
    def _merge_case_collisions(dry_run: bool):
        surviving_tags = []
        for tag_ids in get_case_collision_groups():
            tags = sorted((Tag.objects.get(pk=t) for t in tag_ids),
                          key=lambda t: get_tag_usage(t).total, reverse=True)
            surviving_tag = tags[0]
            for dying_tag in tags[1:]:
                if dry_run:
                    logging.info("Would merge '%s' (%s) into '%s'",
                                 dying_tag, get_tag_usage(dying_tag).description(), surviving_tag)
                else:
                    merge_tag(dying_tag, surviving_tag, admin_bot(), set_nodes_dirty=False)
            if not dry_run:
                surviving_tags.append(surviving_tag)
        if surviving_tags:
            # One bulk pass over every affected node after all the merges, rather than per merge
            logging.info("Setting analysis nodes that use the merged tags dirty")
            set_tag_nodes_dirty(surviving_tags)

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
