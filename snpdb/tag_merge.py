"""
Merging tags together, and the case-collision checks that stop them being created in the first place.

Tag.id is the tag name (CharField primary key) so a merge repoints every foreign key at the surviving tag
then deletes the dying one. Most of those tables have no unique constraint involving tag, so they move in a
single UPDATE - only the two that are unique_together with tag have rows that can't be repointed.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging
from collections import defaultdict
from dataclasses import dataclass, field

from Levenshtein import distance
from django.db import transaction
from django.db.models import Count, Model

from analysis.models import VariantTag
from analysis.models.nodes.filters.tag_node import TagNode, TagNodeTag
from snpdb.models import SampleTag, Tag, TagColor, VCFTag

MERGE_SUGGESTION_MAX_DISTANCE = 1


@dataclass(frozen=True)
class TagUsageEntry:
    """ One foreign key table's use of a tag """
    label: str
    count: int


@dataclass(frozen=True)
class TagUsage:
    tag_id: str
    entries: list[TagUsageEntry]

    @property
    def total(self) -> int:
        return sum(e.count for e in self.entries)

    def description(self) -> str:
        """ Human readable summary of what a merge would move, eg '404 variant tags, 2 analysis nodes' """
        used = [f"{e.count} {e.label}" for e in self.entries if e.count]
        return ", ".join(used) or "nothing"


@dataclass
class TagMergeCounts:
    label: str
    moved: int = 0
    deleted: int = 0


@dataclass
class TagMergeResult:
    dying_tag_id: str
    surviving_tag_id: str
    counts: list[TagMergeCounts] = field(default_factory=list)

    def description(self) -> str:
        parts = []
        for c in self.counts:
            if c.moved or c.deleted:
                parts.append(f"{c.label}: {c.moved} moved, {c.deleted} deleted as duplicates")
        summary = "; ".join(parts) or "nothing was using it"
        return f"Merged '{self.dying_tag_id}' into '{self.surviving_tag_id}' ({summary})"


@dataclass(frozen=True)
class TagForeignKey:
    label: str
    model: type[Model]
    # The other half of a unique_together with tag. Set means rows whose partner already uses the surviving
    # tag can't be repointed, and that the table is small enough to move a row at a time
    unique_with: str = None


# Every FK to Tag
TAG_FOREIGN_KEYS = [
    TagForeignKey("variant tags", VariantTag),
    TagForeignKey("analysis tag nodes", TagNodeTag, unique_with="tag_node_id"),
    TagForeignKey("tag colour settings", TagColor, unique_with="collection_id"),
    TagForeignKey("VCF tags", VCFTag),
    TagForeignKey("sample tags", SampleTag),
]


def get_tag_usage(tag: Tag) -> TagUsage:
    entries = [TagUsageEntry(fk.label, fk.model.objects.filter(tag=tag).count()) for fk in TAG_FOREIGN_KEYS]
    return TagUsage(tag_id=tag.pk, entries=entries)


def get_tag_usage_by_tag_id() -> dict[str, TagUsage]:
    """ Usage for every tag, done as one aggregate per FK table """
    counts_by_model = {}
    for fk in TAG_FOREIGN_KEYS:
        qs = fk.model.objects.values_list("tag_id").annotate(count=Count("pk"))
        counts_by_model[fk.label] = dict(qs.values_list("tag_id", "count"))

    usage = {}
    for tag_id in Tag.objects.values_list("pk", flat=True):
        entries = [TagUsageEntry(fk.label, counts_by_model[fk.label].get(tag_id, 0)) for fk in TAG_FOREIGN_KEYS]
        usage[tag_id] = TagUsage(tag_id=tag_id, entries=entries)
    return usage


def get_case_collision(tag_name: str) -> Tag | None:
    """ Existing tag whose name only differs by case """
    lower = tag_name.lower()
    for tag in Tag.objects.all():
        if tag.pk.lower() == lower and tag.pk != tag_name:
            return tag
    return None


def get_case_collision_groups() -> list[list[str]]:
    """ Groups of >1 tag names that are the same ignoring case """
    by_lower = defaultdict(list)
    for tag_id in Tag.objects.values_list("pk", flat=True):
        by_lower[tag_id.lower()].append(tag_id)
    return [sorted(tag_ids) for tag_ids in by_lower.values() if len(tag_ids) > 1]


def get_merge_suggestions() -> dict[str, list[str]]:
    """ tag_id -> other tags that look like typos of it (case insensitive Levenshtein distance <= 1) """
    tag_ids = list(Tag.objects.values_list("pk", flat=True))
    suggestions = defaultdict(list)
    for i, tag_id in enumerate(tag_ids):
        for other_id in tag_ids[i + 1:]:
            if distance(tag_id.lower(), other_id.lower()) <= MERGE_SUGGESTION_MAX_DISTANCE:
                suggestions[tag_id].append(other_id)
                suggestions[other_id].append(tag_id)
    return {tag_id: sorted(others) for tag_id, others in suggestions.items()}


def _repoint_tag(fk: TagForeignKey, dying_tag: Tag, surviving_tag: Tag) -> TagMergeCounts:
    """ Move fk's rows from dying_tag to surviving_tag """
    counts = TagMergeCounts(label=fk.label)
    dying_qs = fk.model.objects.filter(tag=dying_tag)
    if fk.unique_with is None:
        counts.moved = dying_qs.update(tag=surviving_tag)
        return counts

    # Being unique_together with tag makes a row whose partner already uses the surviving tag say the same
    # thing twice, so it goes. These are small settings tables, so move a row at a time to let the
    # save()/delete() version bumps that other things cache on happen
    taken = set(fk.model.objects.filter(tag=surviving_tag).values_list(fk.unique_with, flat=True))
    for obj in dying_qs:
        if getattr(obj, fk.unique_with) in taken:
            obj.delete()
            counts.deleted += 1
        else:
            obj.tag = surviving_tag
            obj.save()
            counts.moved += 1
    return counts


def _set_tag_nodes_dirty(tag: Tag):
    """ Nodes filtering on either side of a merge return different variants afterwards """
    for node in TagNode.objects.filter(tagnodetag__tag=tag).distinct():
        node.queryset_dirty = True
        node.save()


def merge_tag(dying_tag: Tag, surviving_tag: Tag) -> TagMergeResult:
    """ Repoint everything using dying_tag at surviving_tag, then delete dying_tag. Cannot be undone.
        Repeated variant tags this leaves behind are indistinguishable from ones that were already there
        - @see the variant_tags delete-duplicates management command """
    if dying_tag.pk == surviving_tag.pk:
        raise ValueError("Cannot merge a tag into itself")

    result = TagMergeResult(dying_tag_id=dying_tag.pk, surviving_tag_id=surviving_tag.pk)
    with transaction.atomic():
        for fk in TAG_FOREIGN_KEYS:
            result.counts.append(_repoint_tag(fk, dying_tag, surviving_tag))
        # Everything now points at the surviving tag, so this catches nodes that used either side
        _set_tag_nodes_dirty(surviving_tag)
        dying_tag.delete()

    logging.info(result.description())
    return result
