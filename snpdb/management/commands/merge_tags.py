"""
Merge tags together, and clean up duplicate variant tags.

Tag names are their primary key, so deployments accumulate case variants of the same tag (eg 'artefact' and
'Artefact') that split reporting and node filters. This merges them, and optionally drops variant tag rows
where the same user tagged the same variant with the same tag in the same analysis more than once.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging

from django.core.management import BaseCommand, CommandError

from snpdb.models import Tag
from snpdb.tag_merge import (
    delete_duplicate_variant_tags,
    get_case_collision_groups,
    get_tag_usage,
    merge_tag,
)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--dying', help="Tag to merge away and delete")
        parser.add_argument('--surviving', help="Tag to keep (required with --dying)")
        parser.add_argument('--case-collisions', action='store_true',
                            help="Merge every set of tags that differ only by case, keeping the most used")
        parser.add_argument('--delete-duplicate-variant-tags', action='store_true',
                            help="Delete variant tags repeating the same variant/tag/analysis/user")
        parser.add_argument('--dry-run', action='store_true', help="Report what would change without doing it")

    def handle(self, *args, **options):
        dying_name = options["dying"]
        surviving_name = options["surviving"]
        dry_run = options["dry_run"]

        if bool(dying_name) != bool(surviving_name):
            raise CommandError("--dying and --surviving must be used together")

        if dying_name:
            self._merge(self._get_tag(dying_name), self._get_tag(surviving_name), dry_run)

        if options["case_collisions"]:
            for tag_ids in get_case_collision_groups():
                tags = sorted((Tag.objects.get(pk=t) for t in tag_ids),
                              key=lambda t: get_tag_usage(t).total, reverse=True)
                surviving_tag = tags[0]
                for dying_tag in tags[1:]:
                    self._merge(dying_tag, surviving_tag, dry_run)

        if options["delete_duplicate_variant_tags"]:
            if dry_run:
                logging.info("Would delete duplicate variant tags")
            else:
                num_deleted = delete_duplicate_variant_tags()
                logging.info("Deleted %d duplicate variant tags", num_deleted)

    @staticmethod
    def _get_tag(tag_name: str) -> Tag:
        try:
            return Tag.objects.get(pk=tag_name)
        except Tag.DoesNotExist as dne:
            raise CommandError(f"No such tag: '{tag_name}'") from dne

    @staticmethod
    def _merge(dying_tag: Tag, surviving_tag: Tag, dry_run: bool):
        if dry_run:
            usage = get_tag_usage(dying_tag)
            logging.info("Would merge '%s' (%s) into '%s'", dying_tag, usage.description(), surviving_tag)
        else:
            merge_tag(dying_tag, surviving_tag)
