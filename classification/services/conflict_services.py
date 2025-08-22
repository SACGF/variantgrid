import itertools
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from enum import StrEnum
from functools import cached_property
from typing import Optional, Iterable, Callable, Any
from datetime import datetime
from django.contrib.auth.models import User
from django.db import transaction
from more_itertools.recipes import partition
from classification.enums import AlleleOriginBucket, ConflictSeverity, ShareLevel, ConflictType, TestingContextBucket
from classification.models import AlleleOriginGrouping, EvidenceKeyMap, Conflict, ConflictKey, ConflictHistory, \
    ConflictLab, ClassificationGrouping, ClassificationModification
from genes.hgvs import CHGVS
from library.utils import strip_json, first, sort_and_group_by, RowSpanTable
from snpdb.models import Allele, Lab

"""

@dataclass(frozen=True)
class ClinSigData(ConflictDataRow):
    clin_sig: str
    amp_level: int

    @staticmethod
    def from_data(lab_id: int, share_level: ShareLevel, row: dict) -> 'ClinSigData':
        return ClinSigData(lab_id=lab_id, share_level=share_level, clin_sig=row.get("clinical_significance"), amp_level=row.get("amp_level"))

    def to_json(self) -> dict:
        return {
            "lab_id": self.lab_id,
            "clinical_significance": self.clin_sig,
            "amp_level": self.amp_level
        }
"""

@dataclass
class ConflictDataRow:
    lab_id: int
    share_level: ShareLevel
    classification: Optional[str]
    clinical_significance: Optional[str]
    amp_level: Optional[str]
    exclude: bool = False
    c_hgvs: Optional[CHGVS] = None
    message: Optional[str] = None

    @cached_property
    def lab(self):
        # Lab uses Caching Object Manager
        return Lab.objects.get(id=self.lab_id)

    def can_view(self, user: User) -> bool:
        return self.share_level.has_access(self.lab, user)

    def with_exclude(self, exclude: bool) -> 'ConflictDataRow':
        if self.exclude == exclude:
            return self
        return ConflictDataRow(
            lab_id=self.lab_id,
            share_level=self.share_level,
            classification=self.classification,
            clinical_significance=self.clinical_significance,
            amp_level=self.amp_level,
            exclude=exclude
        )

    def to_json(self) -> dict:
        return strip_json({
            "lab_id": self.lab_id,
            "share_level": self.share_level.value,
            "classification": self.classification,
            "clinical_significance": self.clinical_significance,
            "amp_level": self.amp_level,
            "exclude": self.exclude
        })

    @staticmethod
    def from_json(row: dict) ->'ConflictDataRow' :
        return ConflictDataRow(
            lab_id=row.get("lab_id"),
            share_level=ShareLevel(row.get("share_level")),
            classification=row.get("classification"),
            clinical_significance=row.get("clinical_significance"),
            amp_level=row.get("amp_level"),
            exclude=row.get("exclude")
        )

    @staticmethod
    def from_data(row: dict, lab_id: int, share_level: ShareLevel, c_hgvs: Optional[CHGVS] = None) -> 'ConflictDataRow':
        return ConflictDataRow(
            lab_id=lab_id,
            share_level=share_level,
            classification=row.get("classification"),
            clinical_significance=row.get("clinical_significance"),
            amp_level=row.get("amp_level"),
            c_hgvs=c_hgvs
        )

    @property
    def bucket(self) -> Optional[int]:
        if classification := self.classification:
            return EvidenceKeyMap.clinical_significance_to_bucket().get(classification)
        return None

    def __lt__(self, other: 'ConflictDataRow') -> bool:
        if self.lab_id != other.lab_id:
            return self.lab_id < other.lab_id
        return self.share_level < other.share_level


class ConflictCalculator(ABC):

    def __init__(self, allele: Allele, conflict_key: ConflictKey, conflict_data: list[ConflictDataRow]):
        self.allele = allele
        self.conflict_key = conflict_key
        self.conflict_data = self.exclude_duplicates(conflict_data)

    @property
    def included_conflict_data(self) -> list[ConflictDataRow]:
        return [cd for cd in self.conflict_data if not cd.exclude]

    @abstractmethod
    def calculate_severity(self) -> ConflictSeverity:
        raise NotImplementedError()

    def should_exclude(self, row: ConflictDataRow) -> bool:
        return not row.share_level.is_discordant_level

    def exclude_duplicates(self, data_list: list[ConflictDataRow]) -> list[ConflictDataRow]:
        by_lab_id: dict[int, list[ConflictDataRow]] = defaultdict(list)
        for row in data_list:
            by_lab_id[row.lab_id].append(row)

        updated_list: list[ConflictDataRow] = []
        # TODO use Lab proper ordering
        # but then need to work out efficient way to load Lab object
        for key, values in by_lab_id.items():
            for index, data_row in enumerate(sorted(list(values), reverse=True)):
                if index == 0 and not self.should_exclude(data_row):
                    updated_list.append(data_row)
                else:
                    updated_list.append(data_row.with_exclude(True))
        return list(sorted(updated_list))

    @cached_property
    def involved_lab_ids(self) -> set[int]:
        lab_ids = set()
        for row in self.conflict_data:
            lab_ids.add(row.lab_id)
        return lab_ids

    @cached_property
    def active_labs(self) -> dict[Lab, bool]:
        lab_ids = defaultdict(lambda: False)
        for row in self.conflict_data:
            lab_ids[row.lab] |= not row.exclude
        return dict(lab_ids)

    def data_json(self) -> dict:
        data_rows = []
        included, excluded = partition(lambda x: x.exclude, sorted(self.conflict_data))
        for row in list(included) + list(excluded):
            data_rows.append(row.to_json())
        return {
            "rows": data_rows
        }

    @transaction.atomic
    def log_history(self, override_date: Optional[datetime] = None):
        allele = self.allele
        conflict_key = self.conflict_key
        conflict: Conflict
        created: bool
        if not bool(self.included_conflict_data):
            # if our calculation is empty (e.g. we worked out there's no contributing labs)
            # don't make the conflict just for it to be empty
            # BUT still need to call log_history in case there was data once and now we need to log
            # it's no longer there
            conflict = Conflict.objects.filter(
                allele=allele,
                conflict_type=conflict_key.conflict_type,
                allele_origin_bucket=conflict_key.allele_origin_bucket,
                testing_context_bucket=conflict_key.testing_context_bucket,
                tumor_type_category=conflict_key.tumor_type_category,
            ).first()
            if not conflict:
                return
            created = False
        else:
            conflict, created = Conflict.objects.get_or_create(
                allele=allele,
                conflict_type=conflict_key.conflict_type,
                allele_origin_bucket=conflict_key.allele_origin_bucket,
                testing_context_bucket=conflict_key.testing_context_bucket,
                tumor_type_category=conflict_key.tumor_type_category,
            )

        latest: Optional[ConflictHistory] = None
        if not created:
            latest = ConflictHistory.objects.filter(conflict=conflict, is_latest=True).first()

        # also need to create ConflictLab
        active_labs = self.active_labs
        ConflictLab.objects.bulk_create([
            ConflictLab(
                conflict=conflict,
                lab=lab,
                active=active
            ) for lab, active in active_labs.items()
        ], update_conflicts=True, update_fields=["active"], unique_fields=["conflict", "lab_id"])
        # TODO have something different from withdrawn
        ConflictLab.objects.filter(conflict=conflict).exclude(lab__in=active_labs.keys()).update(active=False)

        c_hgvs_source = self.included_conflict_data
        if not c_hgvs_source:
            c_hgvs_source = self.conflict_data

        c_hgvs_set = [row.c_hgvs for row in c_hgvs_source]
        if c_hgvs_list := [c_hgvs.to_json_short() for c_hgvs in sorted(set(c_hgvs_set)) if c_hgvs]:
            conflict.meta_data["c_hgvs"] = c_hgvs_list
            conflict.save(update_fields=["meta_data"], update_modified=False)

        data = self.data_json()
        if latest:
            if latest.data == data:
                return  # data is already up to date
            else:
                latest.is_latest = False
                latest.save(update_fields=["is_latest"])

        ch = ConflictHistory(
            conflict=conflict,
            data=data,
            severity=self.calculate_severity(),
            is_latest=True
        )
        ch.save()
        if override_date:
            ch.created = override_date
            ch.modified = override_date
            ch.save(update_fields=["created", "modified"], update_modified=False)


class OncPathCalculator(ConflictCalculator):

    def should_exclude(self, row: ConflictDataRow):
        return super().should_exclude(row) or row.bucket is None

    def calculate_severity(self) -> ConflictSeverity:
        oncpath_data = self.included_conflict_data
        if len(oncpath_data) == 0:
            return ConflictSeverity.NO_SUBMISSIONS
        elif len(oncpath_data) == 1:
            return ConflictSeverity.SINGLE_SUBMISSION
        else:
            buckets: set = set()
            values: set = set()
            for row in oncpath_data:
                if bucket := row.bucket:
                    buckets.add(bucket)
                values.add(row.classification)

            if len(buckets) > 1:
                if 3 in buckets:  # don't like this hardcoding, 3 = LP, P
                    return ConflictSeverity.MAJOR
                else:
                    return ConflictSeverity.MEDIUM
            if len(values) > 1:
                return ConflictSeverity.MINOR
            else:
                return ConflictSeverity.SAME


class ClinSigCalculator(ConflictCalculator):

    def should_exclude(self, row: ConflictDataRow):
        return super().should_exclude(row) or row.clinical_significance is None

    def calculate_severity(self) -> ConflictSeverity:
        clin_sig_data = self.included_conflict_data
        if len(clin_sig_data) == 0:
            return ConflictSeverity.NO_SUBMISSIONS
        elif len(clin_sig_data) == 1:
            return ConflictSeverity.SINGLE_SUBMISSION
        else:
            has_tier_1_and_2 = False
            tiers = set()
            for row in self.conflict_data:
                # TODO this is hardcoding values expected in evidence keys
                cs = row.clinical_significance
                if cs == "tier_1_or_2":
                    has_tier_1_and_2 = True
                else:
                    tiers.add(cs)

            # for now ever difference in tier is a MAJOR diff
            # the only exception being tier_1 vs tier_1_or_2 or tier_2 vs tier_1_or_2 which
            # is considered minor
            if len(tiers) > 1:
                return ConflictSeverity.MAJOR
            elif len(tiers) == 1 and has_tier_1_and_2:
                only_tier = first(tiers)
                if only_tier in {"tier_1", "tier_2"}:
                    return ConflictSeverity.MINOR
                else:
                    return ConflictSeverity.MAJOR
            return ConflictSeverity.SAME


def calculate_and_apply_conflicts_for(allele_origin_grouping: AlleleOriginGrouping):
    oncpath_data: list[ConflictDataRow] = []
    clinsig_data: list[ConflictDataRow] = []
    allele_origin_bucket = allele_origin_grouping.allele_origin_bucket

    for cg in allele_origin_grouping.classificationgrouping_set.select_related("latest_classification_modification__classification"):
        lab_id = cg.latest_classification_modification.classification.lab_id
        summary_info = cg.latest_classification_modification.classification.summary_typed
        c_hgvs = cg.latest_classification_modification.imported_c_hgvs_obj

        oncpath_data.append(ConflictDataRow.from_data(row=summary_info.get("pathogenicity", {}), lab_id=lab_id, share_level=cg.share_level_obj, c_hgvs=c_hgvs))
        clinsig_data.append(ConflictDataRow.from_data(row=summary_info.get("somatic", {}), lab_id=lab_id, share_level=cg.share_level_obj, c_hgvs=c_hgvs))

    allele = allele_origin_grouping.allele_grouping.allele
    if allele_origin_bucket == AlleleOriginBucket.SOMATIC:
        conflict_key = ConflictKey(
            conflict_type=ConflictType.CLIN_SIG,
            allele_origin_bucket=allele_origin_bucket,
            testing_context_bucket=allele_origin_grouping.testing_context_bucket,
            tumor_type_category=allele_origin_grouping.tumor_type_category
        )
        ClinSigCalculator(allele, conflict_key, clinsig_data).log_history()

    conflict_key = ConflictKey(
        conflict_type=ConflictType.ONCPATH,
        allele_origin_bucket=allele_origin_bucket,
        testing_context_bucket=allele_origin_grouping.testing_context_bucket,
        tumor_type_category=allele_origin_grouping.tumor_type_category
    )
    OncPathCalculator(allele, conflict_key, oncpath_data).log_history()


def latest_classifications_for_conflict(conflict: Conflict) -> Optional[ClassificationModification]:
    # easy case before we get into
    if grouping := AlleleOriginGrouping.objects.filter(
        allele_grouping__allele=conflict.allele,
        allele_origin_bucket=conflict.allele_origin_bucket,
        testing_context_bucket=conflict.testing_context_bucket,
        tumor_type_category=conflict.tumor_type_category
    ).first():
        # TODO filter out records that don't have a relevant value? Still probably pretty useful to include them
        return [cg.latest_classification_modification for cg in ClassificationGrouping.objects.filter(allele_origin_grouping=grouping)]
    else:
        raise ValueError(f"Can't find AlleleOriginGrouping for conflict {conflict}")


@dataclass(frozen=True)
class ConflictMerge:
    conflicts: list[Conflict]

    @cached_property
    def conflict(self) -> Conflict:
        return first(self.conflicts)

    @cached_property
    def conflict_types(self) -> list[ConflictType]:
        conflict_types = set()
        for conflict in self.conflicts:
            conflict_types.add(ConflictType(conflict.conflict_type))
        return list(sorted(conflict_types))

    @cached_property
    def conflict_type_labels(self) -> list[str]:
        # conflict type label can change depending on allele origin (even if calculated the same)
        return [conflict_type.label_for_context(self.conflict.allele_origin_bucket) for conflict_type in self.conflict_types]

    @staticmethod
    def merge(conflicts: Iterable[Conflict]) -> list['ConflictMerge']:
        # groupby conflict key but with testing context hardcoded to None, so records with the same everything else will be grouped together
        all_merged: list[ConflictMerge] = []
        for group, conflicts in itertools.groupby(sorted(conflicts), lambda c: ConflictKey(
                None,
                c.allele_origin_bucket,
                c.testing_context_bucket,
                c.tumor_type_category,
                c.latest.severity
            )):
            conflict_merge = ConflictMerge(list(conflicts))
            all_merged.append(conflict_merge)
        return all_merged

    def __getattr__(self, method_name):
        # delegate to first conflict object if we get a method we're not familiar with
        def method(*args, **kwargs):
            result = getattr(self.conflict, method_name)
            if isinstance(result, Callable):
                return result(*args, **kwargs)
            else:
                return result

        return method


@staticmethod
def group_conflicts(conflicts: Iterable[Conflict]) -> 'RowSpanTable':
    conflicts = sorted(conflicts, key=lambda c: c.conflict_type)

    table = RowSpanTable(5)

    for allele_origin, sub_conflicts in sort_and_group_by(conflicts, lambda c: c.allele_origin_bucket):
        table.add_cell(0, AlleleOriginBucket(allele_origin).label)

        for testing_context, sub_sub_conflicts in sort_and_group_by(sub_conflicts, lambda c: c.testing_context_bucket):

            table.add_cell(1, TestingContextBucket(testing_context).label)

            for condition, sub_sub_sub_conflicts in sort_and_group_by(sub_sub_conflicts, lambda c: c.tumor_type_category or "Undetermined condition"):

                table.add_cell(2, condition)

                for value_type, sub_sub_sub_sub_conflicts in sort_and_group_by(sub_sub_sub_conflicts, lambda c: ConflictType(c.conflict_type).label):
                    sub_sub_sub_sub_conflicts_list = list(sub_sub_sub_sub_conflicts)
                    if len(sub_sub_sub_sub_conflicts_list) != 1:
                        raise ValueError(f"Multiple conflicts found for {allele_origin} {testing_context} {condition} {sub_sub_sub_sub_conflicts_list}")

                    table.add_cell(3, value_type)
                    table.add_cell(4, sub_sub_sub_sub_conflicts_list[0].latest.get_severity_display)

    return table
