"""
Merging tags together, and the case-collision checks that stop them being created in the first place.

Tag.id is the tag name (CharField primary key) so a merge repoints every foreign key at the surviving tag
then deletes the dying one. Repointing can collide with rows that already use the survivor - those are
duplicates, so we keep the earliest and drop the rest.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field

from Levenshtein import distance
from django.db import transaction
from django.db.models import Count, Model
from django.db.models.signals import post_delete

from analysis.models import VariantTag
from analysis.models.nodes.filters.tag_node import TagNode, TagNodeTag
from analysis.signals.signal_handlers import variant_tag_delete
from snpdb.models import SampleTag, Tag, TagColor, VCFTag

MERGE_SUGGESTION_MAX_DISTANCE = 1
BATCH_SIZE = 1000


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


# Every FK to Tag, with the other fields that make a row a duplicate once it's been repointed
_TAG_FOREIGN_KEYS = [
    ("variant tags", VariantTag, ["variant_id", "analysis_id", "user_id"]),
    ("analysis tag nodes", TagNodeTag, ["tag_node_id"]),
    ("tag colour settings", TagColor, ["collection_id"]),
    ("VCF tags", VCFTag, ["vcf_id"]),
    ("sample tags", SampleTag, ["sample_id"]),
]


def get_tag_usage(tag: Tag) -> TagUsage:
    entries = [TagUsageEntry(label, model.objects.filter(tag=tag).count())
               for label, model, _ in _TAG_FOREIGN_KEYS]
    return TagUsage(tag_id=tag.pk, entries=entries)


def get_tag_usage_by_tag_id() -> dict[str, TagUsage]:
    """ Usage for every tag, done as one aggregate per FK table """
    counts_by_model = {}
    for label, model, _ in _TAG_FOREIGN_KEYS:
        qs = model.objects.values_list("tag_id").annotate(count=Count("pk"))
        counts_by_model[label] = dict(qs.values_list("tag_id", "count"))

    usage = {}
    for tag_id in Tag.objects.values_list("pk", flat=True):
        entries = [TagUsageEntry(label, counts_by_model[label].get(tag_id, 0))
                   for label, _, _ in _TAG_FOREIGN_KEYS]
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


@contextmanager
def _variant_tag_delete_signal_disconnected():
    """ Every VariantTag we delete leaves an identical row behind, so the analyses they belong to see no
        change - skip the per-row 'set tag nodes dirty' + celery task the delete signal would fire.
        Merges bump the affected nodes once at the end instead """
    post_delete.disconnect(variant_tag_delete, sender=VariantTag)
    try:
        yield
    finally:
        post_delete.connect(variant_tag_delete, sender=VariantTag)


def _repoint_tag(model: type[Model], dying_tag: Tag, surviving_tag: Tag,
                 duplicate_fields: list[str], label: str) -> TagMergeCounts:
    """ Move model's rows from dying_tag to surviving_tag, deleting any that would duplicate an existing row.
        Earliest row (lowest pk) wins, so re-tagging history keeps its original event """
    counts = TagMergeCounts(label=label)
    existing_keys = set(model.objects.filter(tag=surviving_tag).values_list(*duplicate_fields))

    move_ids = []
    delete_ids = []
    for row in model.objects.filter(tag=dying_tag).order_by("pk").values_list("pk", *duplicate_fields):
        pk, key = row[0], row[1:]
        if key in existing_keys:
            delete_ids.append(pk)
        else:
            existing_keys.add(key)
            move_ids.append(pk)

    if delete_ids:
        model.objects.filter(pk__in=delete_ids).delete()
        counts.deleted = len(delete_ids)
    if move_ids:
        counts.moved = model.objects.filter(pk__in=move_ids).update(tag=surviving_tag)
    return counts


def _set_tag_nodes_dirty(tag_ids: list[str]):
    """ Nodes filtering on either side of a merge return different variants afterwards """
    nodes_qs = TagNode.objects.filter(tagnodetag__tag__in=tag_ids).distinct()
    for node in nodes_qs:
        node.queryset_dirty = True
        node.save()


def merge_tag(dying_tag: Tag, surviving_tag: Tag) -> TagMergeResult:
    """ Repoint everything using dying_tag at surviving_tag, then delete dying_tag. Cannot be undone """
    if dying_tag.pk == surviving_tag.pk:
        raise ValueError("Cannot merge a tag into itself")

    result = TagMergeResult(dying_tag_id=dying_tag.pk, surviving_tag_id=surviving_tag.pk)
    with transaction.atomic(), _variant_tag_delete_signal_disconnected():
        for label, model, duplicate_fields in _TAG_FOREIGN_KEYS:
            result.counts.append(_repoint_tag(model, dying_tag, surviving_tag, duplicate_fields, label))
        _set_tag_nodes_dirty([surviving_tag.pk])
        dying_tag.delete()

    logging.info(result.description())
    return result


def delete_duplicate_variant_tags(tag: Tag = None) -> int:
    """ Same user tagging the same variant with the same tag in the same analysis more than once.
        Different users re-tagging is agreement data, so is kept """
    qs = VariantTag.objects.all()
    if tag:
        qs = qs.filter(tag=tag)

    seen = set()
    delete_ids = []
    for pk, key in ((row[0], row[1:]) for row in
                    qs.order_by("pk").values_list("pk", "variant_id", "tag_id", "analysis_id", "user_id")):
        if key in seen:
            delete_ids.append(pk)
        else:
            seen.add(key)

    with _variant_tag_delete_signal_disconnected():
        for i in range(0, len(delete_ids), BATCH_SIZE):
            VariantTag.objects.filter(pk__in=delete_ids[i:i + BATCH_SIZE]).delete()
    return len(delete_ids)
