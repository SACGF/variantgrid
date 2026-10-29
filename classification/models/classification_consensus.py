"""
ClassificationConsensus: copying content from previous classifications into a record. Lists the candidates a
curator may copy from, grouped by what copying them would write (ConsensusGroup, compared side by side by
ConsensusComparison), and applies one (ClassificationConsensus.apply_to, which patches as
SubmissionSource.CONSENSUS). Entry points: ClassificationConsensus.allele_consensus_groups,
ClassificationConsensus.gene_consensus_groups and ClassificationConsensus.consensus_patch.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import TYPE_CHECKING, Any, Optional, Union

from classification.enums import CopyScope, SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.models.evidence_key import EvidenceKeyMap, VCPatch
from library.utils import nest_dict
from snpdb.enums import AlleleOriginBucket
from snpdb.user_settings_manager import UserSettingsManager

if TYPE_CHECKING:
    from collections.abc import Iterable

    from django.contrib.auth.models import User

    from classification.models.condition_resolved import ConditionResolved
    from genes.models import GeneSymbol
    from snpdb.models.models_variant import Allele

COPY_SCOPES_ALL = frozenset({CopyScope.ALLELE, CopyScope.GENE})
COPY_SCOPES_GENE = frozenset({CopyScope.GENE})
CONSENSUS_GROUP_LIMIT = 10
EXTERNAL_CANDIDATE_LIMIT = 10

BUCKET_ALLELE_ORIGIN = {  # what the allele_origin evidence key is seeded with per bucket
    AlleleOriginBucket.GERMLINE: "germline",
    AlleleOriginBucket.SOMATIC: "somatic",
}


@dataclass(frozen=True)
class ContentValue:
    """ One copied key's content on a candidate - the thing a row is chosen for. `current` is what the record
        being filled already holds, set only where there is one (the form's Gene Content dialog) """
    key: str
    label: str
    value: Any
    note: Optional[str] = None
    current: Any = None


@dataclass(frozen=True)
class ConsensusGroup:
    """
    Classifications whose copied content is identical - they were copied from each other, so they are one
    choice rather than many. The others are kept so the curator can see the spread of conditions the content
    has been used for, which is usually what decides whether it applies here.
    """
    representative: ClassificationModification
    others: list[ClassificationModification]
    content: list[ContentValue]
    differs_in: list[str] = field(default_factory=list)  # labels of the fields that differ from the first group

    @property
    def other_count(self) -> int:
        return len(self.others)

    @property
    def record_count(self) -> int:
        return len(self.others) + 1

    @property
    def records(self) -> list[ClassificationModification]:
        return [self.representative, *self.others]

    @property
    def content_by_key(self) -> dict[str, ContentValue]:
        return {value.key: value for value in self.content}

    @cached_property
    def conditions(self) -> list[ConditionResolved]:
        """ The distinct conditions this content has been used for, most used first - the tumour type
            spread is usually what decides whether a row applies to the case in hand. A record with no condition
            adds nothing (its fallback is a reference with neither term nor name) """
        counts: Counter[ConditionResolved] = Counter()
        for modification in self.records:
            condition = modification.condition_resolution_obj_fallback
            if condition.plain_text or any(reference.term or reference.name for reference in condition.references):
                counts[condition] += 1
        return [condition for condition, _ in counts.most_common()]


@dataclass(frozen=True)
class ComparisonRow:
    """ One copied key across the groups being compared - a cell per group, None where that group has nothing """
    key: str
    label: str
    cells: list[Optional[ContentValue]]


@dataclass(frozen=True)
class ConsensusComparison:
    """
    The groups side by side, the way the Classify & Report dialog and the create page offer the allele's
    previous classifications: what applying one would write, with the fields the candidates disagree on
    separated from the ones they all carry the same. Nothing about the records themselves (lab id, who
    curated it) is what is being chosen - that is provenance, and the row keeps it to a line.
    """
    groups: list[ConsensusGroup]
    differing: list[ComparisonRow]
    shared: list[ComparisonRow]

    @property
    def copied_count(self) -> int:
        return len(self.differing) + len(self.shared)

    @staticmethod
    def of(groups: list[ConsensusGroup]) -> ConsensusComparison:
        by_key = [group.content_by_key for group in groups]
        keys_in_order: list[str] = []
        labels: dict[str, str] = {}
        for group in groups:
            for value in group.content:
                if value.key not in labels:
                    keys_in_order.append(value.key)
                    labels[value.key] = value.label
        differing = []
        shared = []
        for key in keys_in_order:
            cells = [content.get(key) for content in by_key]
            row = ComparisonRow(key=key, label=labels[key], cells=cells)
            distinct = {(cell.value, cell.note) if cell else None for cell in cells}
            (differing if len(distinct) > 1 else shared).append(row)
        return ConsensusComparison(groups=groups, differing=differing, shared=shared)


@dataclass
class ClassificationConsensus:
    modification: ClassificationModification
    copy_scopes: frozenset[CopyScope] = COPY_SCOPES_ALL

    @staticmethod
    def default_allele_origin_bucket(user: User) -> Optional[AlleleOriginBucket]:
        """ The bucket a new record starts in, from the user's allele origin focus """
        allele_origin_filter = UserSettingsManager.get_user_settings(user).allele_origin_focus
        if allele_origin_filter.value in (AlleleOriginBucket.GERMLINE.value, AlleleOriginBucket.SOMATIC.value):
            return AlleleOriginBucket(allele_origin_filter.value)
        return None

    @staticmethod
    def copy_keys(copy_scopes: frozenset[CopyScope], allele_origin_bucket: Optional[AlleleOriginBucket]) -> list[str]:
        """ The keys that travel when copying content of these scopes into a record of this bucket """
        return [ekey.key for ekey in EvidenceKeyMap.cached().all_keys
                if ekey.effective_copy_scope in copy_scopes
                and ekey.copy_allele_origin_enum.can_copy_to(allele_origin_bucket)]

    @staticmethod
    def gene_scope_keys(allele_origin_bucket: Optional[AlleleOriginBucket]) -> list[str]:
        """ The keys that travel when copying gene level content into a record of this bucket """
        return ClassificationConsensus.copy_keys(COPY_SCOPES_GENE, allele_origin_bucket)

    @staticmethod
    def allele_consensus_groups(allele: Allele, user: User,
                                allele_origin_bucket: Optional[AlleleOriginBucket],
                                limit: int = CONSENSUS_GROUP_LIMIT) -> list[ConsensusGroup]:
        """ The allele's previous classifications that may be copied from, one row per distinct set of copied
            values - the content a pick would write, not the records (@see consensus_groups) """
        qs = ClassificationModification.latest_for_user(user=user, allele=allele, published=True,
                                                       exclude_external_labs=True,
                                                       allele_origin_bucket=allele_origin_bucket) \
            .select_related("classification", "classification__lab", "classification__lab__organization")
        keys = ClassificationConsensus.copy_keys(COPY_SCOPES_ALL, allele_origin_bucket)
        return ClassificationConsensus.consensus_groups(qs, keys, limit=limit)

    @staticmethod
    def external_lab_candidates(allele: Allele, user: User,
                                allele_origin_bucket: Optional[AlleleOriginBucket] = None,
                                limit: int = EXTERNAL_CANDIDATE_LIMIT) -> list[ClassificationModification]:
        """ What other labs have curated on this allele - worth seeing while deciding how to curate, but never
            a copy source: their evidence was assembled under a config and assertion method reviewed elsewhere """
        qs = ClassificationModification.latest_for_user(user=user, allele=allele, published=True,
                                                       allele_origin_bucket=allele_origin_bucket) \
            .filter(classification__lab__external=True) \
            .select_related("classification", "classification__lab", "classification__lab__organization")
        return sorted(qs, key=lambda vcm: vcm.curated_date_check, reverse=True)[:limit]

    @staticmethod
    def has_copyable_content(modification: ClassificationModification, keys: list[str]) -> bool:
        """ Whether copying this record would write anything - a record with none of the keys filled is not a choice """
        return bool(ClassificationConsensus._content(modification, keys))

    @staticmethod
    def _content(modification: ClassificationModification, keys: list[str]) -> tuple:
        """ The copied content of a record, as something two records can be compared on """
        evidence = modification.published_evidence or {}
        content = []
        for key in keys:
            blob = evidence.get(key)
            if blob is None:
                continue
            for part in ('value', 'note'):
                part_value = blob.get(part)
                if part_value is None or (isinstance(part_value, list) and not part_value):
                    continue
                content.append((key, part, json.dumps(part_value, sort_keys=True, default=str)))
        return tuple(content)

    @staticmethod
    def gene_consensus_groups(gene_symbol: Union[str, GeneSymbol], user: User,
                              allele_origin_bucket: Optional[AlleleOriginBucket] = None,
                              exclude_modifications: Optional[Iterable[ClassificationModification]] = None,
                              exclude_classification: Optional[Classification] = None,
                              limit: int = CONSENSUS_GROUP_LIMIT) -> list[ConsensusGroup]:
        """
        The gene level content available to copy in a gene, one row per distinct set of values - most records
        in a gene carry the same content because they were copied from each other.

        AMP tiering and therapy content are gene *and* tumour type, so this only ever lists: which row applies
        to the case in hand is the curator's call.
        """
        gene_q = classification_gene_symbol_filter(gene_symbol)
        if gene_q is None:
            return []

        qs = ClassificationModification.latest_for_user(user=user, published=True, exclude_external_labs=True,
                                                       allele_origin_bucket=allele_origin_bucket) \
            .filter(gene_q).distinct() \
            .select_related("classification", "classification__lab", "classification__lab__organization")
        if exclude_modifications:
            qs = qs.exclude(pk__in=[cm.pk for cm in exclude_modifications])
        if exclude_classification:
            qs = qs.exclude(classification=exclude_classification)

        return ClassificationConsensus.consensus_groups(
            qs, ClassificationConsensus.gene_scope_keys(allele_origin_bucket), limit=limit)

    @staticmethod
    def consensus_groups(modifications: Iterable[ClassificationModification], keys: list[str],
                         limit: int = CONSENSUS_GROUP_LIMIT) -> list[ConsensusGroup]:
        """
        The records grouped by what copying them would write (`keys`, value and note), one row per distinct
        content. The most used content first, then newest; each later row names the fields in which it
        differs from the first (`differs_in`). A record with none of the keys filled is not a choice.
        """
        by_content: dict[tuple, list[ClassificationModification]] = defaultdict(list)
        for cm in modifications:
            if content := ClassificationConsensus._content(cm, keys):
                by_content[content].append(cm)

        e_keys = EvidenceKeyMap.cached()
        groups = []
        for members in by_content.values():
            members.sort(key=lambda vcm: vcm.curated_date_check, reverse=True)
            groups.append(ConsensusGroup(
                representative=members[0], others=members[1:],
                content=ClassificationConsensus._content_values(members[0], keys, e_keys)))
        # The widest used content first - a summary gets re-worded over time, and the many near copies
        # of it are told apart by what differs from the first row, not by reading each one
        groups.sort(key=lambda group: (group.record_count, group.representative.curated_date_check), reverse=True)
        groups = groups[:limit]
        if len(groups) > 1:
            first = {value.key: (value.value, value.note) for value in groups[0].content}
            labels = {value.key: value.label for group in groups for value in group.content}
            for index, group in enumerate(groups[1:], start=1):
                this = {value.key: (value.value, value.note) for value in group.content}
                differing = [key for key in keys if first.get(key) != this.get(key)]
                groups[index] = replace(group, differs_in=[labels[key] for key in differing])
        return groups

    @staticmethod
    def _content_values(modification: ClassificationModification, keys: list[str],
                        e_keys: EvidenceKeyMap) -> list[ContentValue]:
        """ The copied content of a record as a row shows it - labelled, options resolved to their labels """
        evidence = modification.published_evidence or {}
        values = []
        for key in keys:
            blob = evidence.get(key)
            if blob is None:
                continue
            e_key = e_keys.get(key)
            value = e_key.pretty_value(blob)
            note = blob.get('note') or None
            if value in (None, '', []) and not note:
                continue
            values.append(ContentValue(key=key, label=e_key.pretty_label, value=value, note=note))
        return values

    @property
    def allele_origin_bucket(self) -> Optional[AlleleOriginBucket]:
        """ The bucket the copy is going into - a candidate is only ever offered for its own bucket """
        return self.modification.classification.allele_origin_bucket

    def apply_to(self, classification: Classification, user: User):
        """ Copy this record's values into another one, filling the fields it has left empty """
        classification.patch_value(
            patch=self.consensus_patch,
            clear_all_fields=False,
            user=user,
            source=SubmissionSource.CONSENSUS,
            leave_existing_values=True,
            save=True,
            make_patch_fields_immutable=False)
        classification.publish_latest(user)

    @cached_property
    def consensus_patch(self) -> VCPatch:
        evidence = self.modification.published_evidence
        consensus: dict[str, Any] = {}

        # default allele origin - don't use copy consensus because that would copy "likely somatic" etc
        allele_origin_bucket = self.allele_origin_bucket
        if allele_origin := BUCKET_ALLELE_ORIGIN.get(allele_origin_bucket):
            consensus["allele_origin.value"] = allele_origin

        for key in ClassificationConsensus.copy_keys(self.copy_scopes, allele_origin_bucket):
            for part in ['value', 'note']:
                blob = evidence.get(key)
                if blob is None:
                    continue
                part_value = blob.get(part)
                if part_value is None or (isinstance(part_value, list) and not part_value):
                    continue
                consensus[key + '.' + part] = part_value

        patch: dict[str, dict[str, Any]] = nest_dict(consensus)
        return patch
