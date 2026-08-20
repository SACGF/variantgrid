"""
Admin operations on the tag vocabulary - merging, retiring, reinstating, setting allele origin - and the
case-collision checks that stop confusable tags being created in the first place.

Tag.id is the tag name (CharField primary key) so a merge repoints every foreign key at the surviving tag
then retires the dying one. Variant tags have no unique constraint involving tag so they move in a single
UPDATE - the two settings tables that are unique_together with tag have rows that can't be repointed.

Retiring rather than deleting means an old tag name still resolves, so "where did my tag go?" has an answer.
Each operation writes an auditlog LogEntry. Tag isn't registered with auditlog - we only want the deliberate
operations, not every save - so the entries are written by hand, and the bulk repointing counts have to be
passed in because a queryset .update() fires no signals for auditlog to see.

@see https://github.com/SACGF/variantgrid/issues/1751
"""
import logging
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field

from Levenshtein import distance
from auditlog.models import LogEntry
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.db import transaction
from django.db.models import Count, F, Model, QuerySet, TextChoices
from django.utils import timezone

from analysis.models import Analysis, NodeStatus, VariantTag
from analysis.models.nodes.analysis_node import AnalysisEdge, AnalysisNode, NodeVersion
from analysis.models.nodes.filters.tag_node import TagNode, TagNodeTag
from classification.enums import AlleleOriginBucket
from snpdb.models import TAG_ALLELE_ORIGIN_CHOICES, Tag, TagColor

MERGE_SUGGESTION_MAX_DISTANCE = 1


class TagOperation(TextChoices):
    """ Stored in LogEntry.additional_data["operation"] - not a model field """
    MERGE = "merge", "Merged"
    RETIRE = "retire", "Retired"
    REINSTATE = "reinstate", "Reinstated"
    SET_ALLELE_ORIGIN = "allele_origin", "Allele origin"


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
    for tag_id in Tag.live_qs().values_list("pk", flat=True):
        entries = [TagUsageEntry(fk.label, counts_by_model[fk.label].get(tag_id, 0)) for fk in TAG_FOREIGN_KEYS]
        usage[tag_id] = TagUsage(tag_id=tag_id, entries=entries)
    return usage


def get_case_collision(tag_name: str) -> Tag | None:
    """ Existing tag whose name only differs by case. Includes retired tags - the name is still taken,
        and the caller wants to tell people where it went """
    lower = tag_name.lower()
    for tag in Tag.objects.all():
        if tag.pk.lower() == lower and tag.pk != tag_name:
            return tag
    return None


def get_case_collision_groups() -> list[list[str]]:
    """ Groups of >1 live tag names that are the same ignoring case """
    by_lower = defaultdict(list)
    for tag_id in Tag.live_qs().values_list("pk", flat=True):
        by_lower[tag_id.lower()].append(tag_id)
    return [sorted(tag_ids) for tag_ids in by_lower.values() if len(tag_ids) > 1]


def get_merge_suggestions() -> dict[str, list[str]]:
    """ tag_id -> other tags that look like typos of it (case insensitive Levenshtein distance <= 1) """
    tag_ids = list(Tag.live_qs().values_list("pk", flat=True))
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


def set_tag_nodes_dirty(tags: Iterable[Tag]):
    """ Nodes filtering on any of these tags return different variants after a merge/retire/reinstate,
        and everything downstream of them re-runs too. Saving a node cascades version-bump saves through
        its descendants a row at a time, so this bumps everything in bulk the way reload_analysis_nodes
        does - the periodic scheduler picks the dirty nodes up from there """
    tag_nodes = list(TagNode.objects.filter(tagnodetag__tag__in=tags).distinct())
    if not tag_nodes:
        return

    analysis_ids = {node.analysis_id for node in tag_nodes}
    children_by_parent = defaultdict(list)
    edges_qs = AnalysisEdge.objects.filter(parent__analysis_id__in=analysis_ids)
    for parent_id, child_id in edges_qs.values_list("parent_id", "child_id"):
        children_by_parent[parent_id].append(child_id)

    dirty_node_ids = {node.pk for node in tag_nodes}
    queue = list(dirty_node_ids)
    while queue:
        for child_id in children_by_parent.get(queue.pop(), []):
            if child_id not in dirty_node_ids:
                dirty_node_ids.add(child_id)
                queue.append(child_id)

    nodes_qs = AnalysisNode.objects.filter(pk__in=dirty_node_ids)
    nodes_qs.update(version=F("version") + 1, status=NodeStatus.DIRTY,
                    count=None, errors=None, cloned_from=None)
    NodeVersion.objects.bulk_create([NodeVersion(node_id=node_id, version=version)
                                     for node_id, version in nodes_qs.values_list("pk", "version")],
                                    ignore_conflicts=True)
    Analysis.objects.filter(pk__in=analysis_ids).update(modified=timezone.now())

    # Tag names appear in auto generated node names, so refresh the tag nodes' labels
    renamed = []
    for node in tag_nodes:
        if node.auto_node_name and (name := node.get_node_name()) != node.name:
            node.name = name
            renamed.append(node)
    if renamed:
        TagNode.objects.bulk_update(renamed, ["name"])


def log_tag_operation(tag: Tag, operation: TagOperation, user: User, **details) -> LogEntry:
    """ Append-only record of a deliberate change to the vocabulary. Creating a tag isn't logged - the Tag
        row is that record - and repointing is done in bulk, so what moved has to be passed in as details """
    return LogEntry.objects.log_create(
        tag,
        force_log=True,
        action=LogEntry.Action.UPDATE,
        actor=user,
        additional_data={"operation": operation.value, **details},
    )


def _allele_origin_label(bucket: str) -> str:
    return dict(TAG_ALLELE_ORIGIN_CHOICES).get(bucket, bucket or "")


def describe_tag_operation(log_entry: LogEntry) -> str:
    """ Human readable summary of one operation, from the details stashed at the time """
    data = log_entry.additional_data or {}
    operation = data.get("operation")
    if operation == TagOperation.MERGE:
        moves = [f"{c['label']}: {c['moved']} moved, {c['deleted']} deleted as duplicates"
                 for c in data.get("counts", []) if c["moved"] or c["deleted"]]
        summary = "; ".join(moves) or "nothing was using it"
        return f"into '{data.get('into')}' ({summary})"
    if operation == TagOperation.RETIRE:
        parts = []
        if reason := data.get("reason"):
            parts.append(reason)
        if still_used_by := data.get("still_used_by"):
            parts.append(f"still used by {still_used_by}")
        return " - ".join(parts)
    if operation == TagOperation.REINSTATE:
        if was_merged_into := data.get("was_merged_into"):
            return f"was merged into '{was_merged_into}' - those rows stay where they went"
        return ""
    if operation == TagOperation.SET_ALLELE_ORIGIN:
        return f"{_allele_origin_label(data.get('from'))} -> {_allele_origin_label(data.get('to'))}"
    return ""


def get_tag_operations(tag: Tag = None) -> QuerySet[LogEntry]:
    """ Operation history, newest first - for one tag, or the whole vocabulary. LogEntry is shared with
        everything auditlog tracks, so this only returns entries we wrote """
    if tag is not None:
        qs = LogEntry.objects.get_for_object(tag)
    else:
        qs = LogEntry.objects.filter(content_type=ContentType.objects.get_for_model(Tag))
    return qs.filter(additional_data__has_key="operation").order_by("-timestamp")


def merge_tag(dying_tag: Tag, surviving_tag: Tag, user: User, set_nodes_dirty: bool = True) -> TagMergeResult:
    """ Repoint everything using dying_tag at surviving_tag, then retire dying_tag. The repointing cannot
        be undone - reinstating the tag afterwards gets the name back, not the rows.
        Repeated variant tags this leaves behind are indistinguishable from ones that were already there
        - @see the variant_tags delete-duplicates management command.
        A caller doing several merges passes set_nodes_dirty=False and calls set_tag_nodes_dirty() itself
        with all the surviving tags at the end, so each affected node only re-runs once """
    if dying_tag.pk == surviving_tag.pk:
        raise ValueError("Cannot merge a tag into itself")
    if not surviving_tag.active:
        raise ValueError(f"Cannot merge into '{surviving_tag}', which is retired")

    result = TagMergeResult(dying_tag_id=dying_tag.pk, surviving_tag_id=surviving_tag.pk)
    with transaction.atomic():
        for fk in TAG_FOREIGN_KEYS:
            result.counts.append(_repoint_tag(fk, dying_tag, surviving_tag))
        if set_nodes_dirty:
            # Everything now points at the surviving tag, so this catches nodes that used either side
            set_tag_nodes_dirty([surviving_tag])
        dying_tag.merged_into = surviving_tag
        dying_tag.retired = timezone.now()
        dying_tag.save()
        log_tag_operation(dying_tag, TagOperation.MERGE, user,
                          into=surviving_tag.pk,
                          counts=[asdict(c) for c in result.counts])

    logging.info(result.description())
    return result


def retire_tag(tag: Tag, user: User, reason: str = None) -> Tag:
    """ Stop offering a tag without touching what already uses it - historical rows and stats keep working.
        Use merge instead when the taggings should move somewhere """
    if not tag.active:
        raise ValueError(f"Tag '{tag}' is already retired")

    usage = get_tag_usage(tag)
    with transaction.atomic():
        tag.retired = timezone.now()
        tag.save()
        # Nodes offering this tag return the same variants, but the tag is no longer selectable
        set_tag_nodes_dirty([tag])
        log_tag_operation(tag, TagOperation.RETIRE, user,
                          reason=reason,
                          still_used_by=usage.description())

    logging.info("Retired tag '%s' (%s)", tag, usage.description())
    return tag


def reinstate_tag(tag: Tag, user: User) -> Tag:
    """ Put a retired tag back into circulation. A merged tag keeps merged_into cleared on the way back so
        it stops claiming its rows live somewhere else - they aren't coming back """
    if tag.active:
        raise ValueError(f"Tag '{tag}' is not retired")

    with transaction.atomic():
        merged_into_id = tag.merged_into_id
        tag.retired = None
        tag.merged_into = None
        tag.save()
        set_tag_nodes_dirty([tag])
        log_tag_operation(tag, TagOperation.REINSTATE, user, was_merged_into=merged_into_id)

    logging.info("Reinstated tag '%s'", tag)
    return tag


def set_tag_allele_origin(tag: Tag, allele_origin_bucket: str, user: User) -> Tag:
    """ Which side of the house a tag belongs to - it's what the tag stats page filters on. Tags stay
        offered everywhere regardless, so nothing already tagged or filtering on the tag changes """
    bucket = AlleleOriginBucket(allele_origin_bucket)
    previous = tag.allele_origin_bucket
    if previous == bucket:
        return tag

    with transaction.atomic():
        tag.allele_origin_bucket = bucket
        tag.save()
        log_tag_operation(tag, TagOperation.SET_ALLELE_ORIGIN, user, **{"from": previous, "to": bucket.value})

    logging.info("Tag '%s' allele origin %s -> %s", tag, previous, bucket.value)
    return tag
