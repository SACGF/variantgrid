"""
Which evidence keys moved when a record was reclassified, towards benign or towards pathogenic.
EvidenceMovements is the entry point (analytics.py:ReclassificationAnalytics.evidence builds it over the filtered
reclassifications); count_evidence_movements does the diff in SQL and collapse_evidence folds the quiet keys
into one row for the chart.
"""
from collections import defaultdict
from functools import cached_property
from typing import Optional

from django.core.cache import cache
from django.db import connection
from django.db.models import QuerySet

from classification.enums.classification_enums import CriteriaEvaluation, SpecialEKeys
from classification.models.classification import ClassificationModification
from classification.models.classification_reclassification_models import ReclassificationEvent
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from classification.services.reclassification_analytics.results import EvidenceMovement

EVIDENCE_MOVEMENT_CACHE_SECONDS = 60 * 60
# the significance itself changes by definition, the dates get stamped alongside it, and the report and
# internal notes keys track a record's paperwork rather than the evidence behind the call
EVIDENCE_KEYS_EXCLUDED = {
    SpecialEKeys.CLINICAL_SIGNIFICANCE,
    SpecialEKeys.CURATION_DATE,
    SpecialEKeys.CURATION_VERIFIED_DATE,
    "internal_use",
    "report_date",
    "report_id",
    "report_type",
}
# criteria always stand on their own, the busiest few other keys join them, and the tail folds into one row
TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT = 8
OTHER_EVIDENCE_GROUP = "__other__"

EvidenceMovementCount = tuple[str, bool, Optional[str], Optional[str], int]
""" (key, towards benign, old criteria value, new criteria value, number of reclassifications) """


class EvidenceMovements:

    def __init__(self, reclassifications_qs: QuerySet[ReclassificationEvent], fold_lab_namespaces: bool,
                 cache_scope: str):
        """
        :param fold_lab_namespaces: a cross lab view counts each lab's version of a criterion as one
        :param cache_scope: identifies the filters, so the cached diff is only reused for the same page
        """
        self.reclassifications_qs = reclassifications_qs
        self.fold_lab_namespaces = fold_lab_namespaces
        self.cache_scope = cache_scope

    @cached_property
    def evidence_key_order(self) -> list[EvidenceKey]:
        """ Criteria lead in ACMG strength order, then the rest as they appear on a classification """
        evidence_keys = EvidenceKeyMap.cached()
        criteria = sorted(evidence_keys.criteria(),
                          key=lambda e_key: (-CriteriaEvaluation.POINTS.get(e_key.default_crit_evaluation, 0),
                                             e_key.pretty_label.lower()))
        criteria_keys = {e_key.key for e_key in criteria}
        return criteria + [e_key for e_key in evidence_keys.all_keys
                           if e_key.key not in criteria_keys
                           and e_key.key not in EVIDENCE_KEYS_EXCLUDED]

    @cached_property
    def criteria_keys(self) -> set[str]:
        return {e_key.key for e_key in EvidenceKeyMap.cached().criteria()}

    @cached_property
    def evidence_group_by_key(self) -> dict[str, str]:
        """
        Criteria are namespaced per lab config, so a cross lab view folds each lab's version of a criterion
        together while a single lab view keeps that lab's own keys.
        """
        if not self.fold_lab_namespaces:
            return {}
        return {e_key.key: e_key.without_namespace for e_key in EvidenceKeyMap.cached().criteria()}

    def evidence_group(self, key: str) -> str:
        return self.evidence_group_by_key.get(key, key)

    @cached_property
    def evidence_group_order(self) -> list[tuple[str, str]]:
        """ (group, label) in display order, labelled by the first key that lands in each group """
        groups = []
        seen = set()
        for e_key in self.evidence_key_order:
            group = self.evidence_group(e_key.key)
            if group not in seen:
                seen.add(group)
                groups.append((group, e_key.pretty_label))
        return groups

    @cached_property
    def towards_pathogenic(self) -> list[EvidenceMovement]:
        return self._movements(towards_benign=False)

    @cached_property
    def towards_benign(self) -> list[EvidenceMovement]:
        return self._movements(towards_benign=True)

    def _movements(self, towards_benign: bool) -> list[EvidenceMovement]:
        applied = defaultdict(int)
        unapplied = defaultdict(int)
        strengthened = defaultdict(int)
        weakened = defaultdict(int)
        changed = defaultdict(int)

        for key, row_towards_benign, old_value, new_value, total in self._movement_counts:
            if row_towards_benign != towards_benign:
                continue
            group = self.evidence_group(key)
            if key not in self.criteria_keys:
                changed[group] += total
                continue

            old_met = CriteriaEvaluation.is_met(old_value)
            new_met = CriteriaEvaluation.is_met(new_value)
            if old_met != new_met:
                (applied if new_met else unapplied)[group] += total
            else:
                old_strength = abs(CriteriaEvaluation.POINTS.get(old_value, 0))
                new_strength = abs(CriteriaEvaluation.POINTS.get(new_value, 0))
                if new_strength > old_strength:
                    strengthened[group] += total
                elif new_strength < old_strength:
                    weakened[group] += total
                else:
                    changed[group] += total

        moved = set(applied) | set(unapplied) | set(strengthened) | set(weakened) | set(changed)
        criteria_groups = {self.evidence_group(key) for key in self.criteria_keys}
        return [EvidenceMovement(key=group, label=label, applied=applied[group], unapplied=unapplied[group],
                                 strengthened=strengthened[group], weakened=weakened[group],
                                 changed=changed[group], is_criteria=group in criteria_groups)
                for group, label in self.evidence_group_order if group in moved]

    @cached_property
    def _movement_counts(self) -> list[EvidenceMovementCount]:
        event_ids = list(self.reclassifications_qs.filter(from_modification__isnull=False)
                         .values_list('pk', flat=True))
        if not event_ids:
            return []

        cache_key = f"reclassification_evidence_movements:{self.cache_scope}:{len(event_ids)}"
        counts = cache.get(cache_key)
        if counts is None:
            counts = count_evidence_movements(event_ids, sorted(self.criteria_keys))
            cache.set(cache_key, counts, EVIDENCE_MOVEMENT_CACHE_SECONDS)
        return counts


def collapse_evidence(movements: list[EvidenceMovement]) -> tuple[list[EvidenceMovement], int]:
    """ Criteria and the busiest few other keys keep their own row, the long tail folds into one """
    others = sorted((movement for movement in movements if not movement.is_criteria),
                    key=lambda movement: -movement.total)
    kept = {movement.key for movement in others[:TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT]}
    folded = others[TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT:]
    rows = [movement for movement in movements if movement.is_criteria or movement.key in kept]
    if folded:
        rows.append(EvidenceMovement(
            key=OTHER_EVIDENCE_GROUP,
            label=f"Other ({len(folded)} keys)",
            applied=sum(movement.applied for movement in folded),
            unapplied=sum(movement.unapplied for movement in folded),
            strengthened=sum(movement.strengthened for movement in folded),
            weakened=sum(movement.weakened for movement in folded),
            changed=sum(movement.changed for movement in folded),
            is_criteria=False))
    return rows, len(folded)


def count_evidence_movements(event_ids: list[int], criteria_keys: list[str]) -> list[EvidenceMovementCount]:
    """
    Diffs each event's two published_evidence blobs in the database - shipping thousands of whole
    evidence documents back to python to compare them would be far more expensive. Criteria carry the
    strength on each side so the caller can tell an on/off from a strength shift, everything else only
    needs to say that it moved.
    """
    event_table = ReclassificationEvent._meta.db_table
    modification_table = ClassificationModification._meta.db_table
    sql = f"""
        SELECT movement.key,
               movement.towards_benign,
               movement.old_value,
               movement.new_value,
               COUNT(*) AS total
        FROM (
            SELECT event.id,
                   entry.key,
                   event.significance_delta > 0 AS towards_benign,
                   CASE WHEN entry.key = ANY(%s) THEN
                       COALESCE(old_version.published_evidence -> entry.key -> 'value',
                                old_version.published_evidence -> entry.key) #>> '{{}}' END AS old_value,
                   CASE WHEN entry.key = ANY(%s) THEN
                       COALESCE(new_version.published_evidence -> entry.key -> 'value',
                                new_version.published_evidence -> entry.key) #>> '{{}}' END AS new_value
            FROM {event_table} event
            JOIN {modification_table} old_version ON old_version.id = event.from_modification_id
            JOIN {modification_table} new_version ON new_version.id = event.to_modification_id
            CROSS JOIN LATERAL (
                SELECT key FROM jsonb_object_keys(COALESCE(old_version.published_evidence, '{{}}'::jsonb)) AS key
                UNION
                SELECT key FROM jsonb_object_keys(COALESCE(new_version.published_evidence, '{{}}'::jsonb)) AS key
            ) entry
            WHERE event.id = ANY(%s)
              AND event.significance_delta <> 0
              AND COALESCE(old_version.published_evidence -> entry.key -> 'value',
                           old_version.published_evidence -> entry.key)
                  IS DISTINCT FROM
                  COALESCE(new_version.published_evidence -> entry.key -> 'value',
                           new_version.published_evidence -> entry.key)
        ) movement
        GROUP BY 1, 2, 3, 4
    """
    with connection.cursor() as cursor:
        cursor.execute(sql, [criteria_keys, criteria_keys, event_ids])
        return [(key, towards_benign, old_value, new_value, total)
                for key, towards_benign, old_value, new_value, total in cursor.fetchall()]
